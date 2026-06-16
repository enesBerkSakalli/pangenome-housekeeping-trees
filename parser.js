/**
 * parser.js — Newick tree parser and layout engine
 * Converts a Newick string to a tree object and computes (x,y) coordinates.
 */

/**
 * Parse a Newick string into a nested node object.
 * Each node: { name, length, children[] }
 */
function parseNewick(s) {
  const tokens = s.match(/([^:;,()\s]+|:[\d.]+|[;,()\s])/g);
  let pos = 0;

  function parseNode() {
    const node = { name: '', length: 0, children: [] };

    if (tokens[pos] === '(') {
      pos++; // consume '('
      node.children.push(parseNode());
      while (tokens[pos] === ',') {
        pos++;
        node.children.push(parseNode());
      }
      pos++; // consume ')'
    }

    // Node name (possibly empty)
    if (pos < tokens.length && tokens[pos] !== ':' && tokens[pos] !== ')' &&
        tokens[pos] !== ',' && tokens[pos] !== ';') {
      node.name = tokens[pos++];
    }

    // Branch length
    if (pos < tokens.length && tokens[pos] && tokens[pos].startsWith(':')) {
      node.length = parseFloat(tokens[pos].slice(1));
      pos++;
    }

    return node;
  }

  return parseNode();
}

/**
 * Compute maximum depth from root to leaf (sum of branch lengths).
 */
function computeMaxDepth(node, depth = 0) {
  if (!node.children.length) return depth + node.length;
  return Math.max(...node.children.map(c => computeMaxDepth(c, depth + node.length)));
}

/**
 * Collect all leaves in left-to-right order.
 */
function getLeaves(node, leaves = []) {
  if (!node.children.length) { leaves.push(node); return leaves; }
  node.children.forEach(c => getLeaves(c, leaves));
  return leaves;
}

/**
 * Assign (x, y) coordinates for RECTANGULAR layout.
 *   x = cumulative branch length from root (horizontal axis)
 *   y = leaf index (vertical axis), internal = average of children
 */
function layoutRectangular(root) {
  const leaves = getLeaves(root);
  const leafCount = leaves.length;
  const maxDepth = computeMaxDepth(root, 0);

  // Assign y positions to leaves
  leaves.forEach((leaf, i) => { leaf.yi = i; });

  function assignCoords(node, depth) {
    node.x = depth + node.length;
    if (!node.children.length) {
      node.y = node.yi;
    } else {
      node.children.forEach(c => assignCoords(c, node.x));
      node.y = (node.children[0].y + node.children[node.children.length - 1].y) / 2;
    }
    // Normalize x to [0,1]
    node.xNorm = maxDepth > 0 ? node.x / maxDepth : 0;
  }

  root.length = 0;
  assignCoords(root, 0);

  // Normalize y to [0,1]
  function normalizeY(node) {
    node.yNorm = leafCount > 1 ? node.y / (leafCount - 1) : 0.5;
    node.children.forEach(normalizeY);
  }
  normalizeY(root);

  return { root, leafCount, maxDepth };
}

/**
 * Assign (angle, radius) coordinates for RADIAL layout.
 */
function layoutRadial(root) {
  const leaves = getLeaves(root);
  const leafCount = leaves.length;
  const maxDepth = computeMaxDepth(root, 0);

  leaves.forEach((leaf, i) => { leaf.yi = i; });

  function assignDepth(node, depth) {
    node.x = depth + node.length;
    node.xNorm = maxDepth > 0 ? node.x / maxDepth : 0;
    if (!node.children.length) {
      node.angle = (node.yi / (leafCount)) * 2 * Math.PI;
    } else {
      node.children.forEach(c => assignDepth(c, node.x));
      const angles = node.children.map(c => c.angle);
      node.angle = (Math.min(...angles) + Math.max(...angles)) / 2;
    }
  }
  root.length = 0;
  assignDepth(root, 0);

  return { root, leafCount, maxDepth };
}

/**
 * Flatten tree to list of nodes for easy iteration.
 */
function flattenTree(node, list = []) {
  list.push(node);
  node.children.forEach(c => flattenTree(c, list));
  return list;
}
