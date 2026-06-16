(function () {
  const data = window.PANGENOME_VARIANT_TREES;
  const genes = Object.keys(data || {});
  if (!genes.length) {
    throw new Error("No 3D variant tree data found.");
  }

  const palette = {
    "Protein-altering": "#ff7a59",
    "Regulatory/splice": "#ffd166",
    "Background": "#4cc9f0",
    "Internal": "#94a3b8",
  };
  const layerColors = ["rgba(20, 50, 79, 0.24)", "rgba(29, 77, 120, 0.20)", "rgba(43, 107, 158, 0.18)"];
  const layerStroke = ["rgba(125, 211, 252, 0.12)", "rgba(125, 211, 252, 0.16)", "rgba(125, 211, 252, 0.20)"];

  const stageEl = document.getElementById("stage");
  const tooltipEl = document.getElementById("tooltip");
  const geneButtonsEl = document.getElementById("gene-buttons");
  const geneMetaEl = document.getElementById("gene-meta");
  const sceneTitleEl = document.getElementById("scene-title");

  const canvas = document.createElement("canvas");
  canvas.className = "stage-canvas";
  stageEl.appendChild(canvas);
  const ctx = canvas.getContext("2d");

  const state = {
    gene: genes[0],
    yaw: -0.7,
    pitch: 0.6,
    zoom: 1.05,
    dragging: false,
    lastX: 0,
    lastY: 0,
    hoverLeaf: null,
  };

  let currentPayload = null;
  let currentNodes = [];
  let currentEdges = [];
  let projectedLeaves = [];

  function renderGene(gene) {
    currentPayload = data[gene];
    currentNodes = Object.values(currentPayload.nodes);
    currentEdges = currentPayload.edges;
    sceneTitleEl.textContent = `${gene} · ${currentPayload.stats.selected_taxa} taxa`;

    const counts = currentPayload.stats.consequence_groups;
    geneMetaEl.innerHTML = [
      row("Description", currentPayload.description),
      row("Taxa", String(currentPayload.stats.selected_taxa)),
      row("Source variants", String(currentPayload.stats.total_variants_in_source)),
      row("Nodes", String(currentPayload.stats.node_count)),
      row("Protein-altering", String(counts["Protein-altering"] || 0)),
      row("Regulatory/splice", String(counts["Regulatory/splice"] || 0)),
      row("Background", String(counts["Background"] || 0)),
    ].join("");

    [...geneButtonsEl.querySelectorAll(".gene-button")].forEach((button) => {
      button.classList.toggle("active", button.dataset.gene === gene);
    });
  }

  function row(label, value) {
    return `<div class="meta-row"><span>${label}</span><strong>${value}</strong></div>`;
  }

  function resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const width = stageEl.clientWidth;
    const height = stageEl.clientHeight;
    canvas.width = Math.round(width * dpr);
    canvas.height = Math.round(height * dpr);
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  function project(point) {
    const yawCos = Math.cos(state.yaw);
    const yawSin = Math.sin(state.yaw);
    const pitchCos = Math.cos(state.pitch);
    const pitchSin = Math.sin(state.pitch);

    let x = point.x;
    let y = -point.y;
    let z = point.z;

    const x1 = x * yawCos - z * yawSin;
    const z1 = x * yawSin + z * yawCos;
    const y1 = y * pitchCos - z1 * pitchSin;
    const z2 = y * pitchSin + z1 * pitchCos;

    const perspective = 680 / (680 - z2);
    const scale = 0.74 * state.zoom * perspective;
    return {
      x: (canvas.clientWidth / 2) + (x1 * scale),
      y: (canvas.clientHeight / 2) + (y1 * scale),
      depth: z2,
      scale,
    };
  }

  function drawPlane(z, index) {
    const corners = [
      { x: -500, y: -300, z },
      { x: 500, y: -300, z },
      { x: 500, y: 300, z },
      { x: -500, y: 300, z },
    ].map(project);

    ctx.beginPath();
    ctx.moveTo(corners[0].x, corners[0].y);
    corners.slice(1).forEach((corner) => ctx.lineTo(corner.x, corner.y));
    ctx.closePath();
    ctx.fillStyle = layerColors[index];
    ctx.strokeStyle = layerStroke[index];
    ctx.lineWidth = 1.2;
    ctx.fill();
    ctx.stroke();
  }

  function updateHover(clientX, clientY) {
    const rect = canvas.getBoundingClientRect();
    const localX = clientX - rect.left;
    const localY = clientY - rect.top;
    let best = null;

    projectedLeaves.forEach((leaf) => {
      const dx = leaf.screen.x - localX;
      const dy = leaf.screen.y - localY;
      const distance = Math.sqrt((dx * dx) + (dy * dy));
      if (distance <= leaf.hitRadius && (!best || distance < best.distance)) {
        best = { leaf, distance };
      }
    });

    state.hoverLeaf = best ? best.leaf : null;
    if (!state.hoverLeaf) {
      tooltipEl.hidden = true;
      return;
    }

    const node = state.hoverLeaf.node;
    tooltipEl.hidden = false;
    tooltipEl.style.left = `${clientX + 14}px`;
    tooltipEl.style.top = `${clientY + 14}px`;
    tooltipEl.innerHTML = [
      `<strong><code>${node.variant_id}</code></strong>`,
      node.rsid ? `<div>${node.rsid}</div>` : "",
      `<div>${node.consequence_group}</div>`,
      `<div>${node.consequence}</div>`,
      `<div>Active populations: ${node.active_populations}</div>`,
      `<div>Layer: ${node.layer_label}</div>`,
    ].join("");
  }

  function animate() {
    if (!currentPayload) {
      requestAnimationFrame(animate);
      return;
    }

    ctx.clearRect(0, 0, canvas.clientWidth, canvas.clientHeight);

    drawPlane(-180, 0);
    drawPlane(0, 1);
    drawPlane(180, 2);

    const projectedById = {};
    projectedLeaves = [];
    currentNodes.forEach((node) => {
      const screen = project(node);
      projectedById[node.id] = screen;
      if (node.is_leaf) {
        projectedLeaves.push({
          node,
          screen,
          hitRadius: Math.max(6, 4.5 * screen.scale),
        });
      }
    });

    currentEdges.forEach((edge) => {
      const source = projectedById[edge.source];
      const target = projectedById[edge.target];
      if (!source || !target) return;
      const alpha = 0.16 + (Math.max(source.depth, target.depth) + 320) / 1200;
      ctx.beginPath();
      ctx.moveTo(source.x, source.y);
      ctx.lineTo(target.x, target.y);
      ctx.strokeStyle = `rgba(143, 168, 200, ${Math.min(0.55, alpha)})`;
      ctx.lineWidth = 1;
      ctx.stroke();
    });

    const nodes = currentNodes
      .map((node) => ({ node, screen: projectedById[node.id] }))
      .sort((a, b) => a.screen.depth - b.screen.depth);

    nodes.forEach(({ node, screen }) => {
      const radius = node.is_leaf
        ? Math.max(1.8, 2.8 * screen.scale)
        : Math.max(1.3, (1.4 + Math.min(3, Math.log2(node.taxon_count + 1))) * screen.scale);
      ctx.beginPath();
      ctx.arc(screen.x, screen.y, radius, 0, Math.PI * 2);
      ctx.fillStyle = node.is_leaf ? palette[node.consequence_group] : palette.Internal;
      ctx.shadowBlur = node.is_leaf ? 14 : 6;
      ctx.shadowColor = ctx.fillStyle;
      ctx.fill();
      ctx.shadowBlur = 0;
    });

    requestAnimationFrame(animate);
  }

  function onPointerDown(event) {
    state.dragging = true;
    state.lastX = event.clientX;
    state.lastY = event.clientY;
  }

  function onPointerMove(event) {
    if (state.dragging) {
      const dx = event.clientX - state.lastX;
      const dy = event.clientY - state.lastY;
      state.yaw += dx * 0.006;
      state.pitch = Math.max(-1.2, Math.min(1.2, state.pitch + (dy * 0.006)));
      state.lastX = event.clientX;
      state.lastY = event.clientY;
      tooltipEl.hidden = true;
      return;
    }
    updateHover(event.clientX, event.clientY);
  }

  function onPointerUp() {
    state.dragging = false;
  }

  function onWheel(event) {
    event.preventDefault();
    state.zoom = Math.max(0.45, Math.min(2.3, state.zoom - (event.deltaY * 0.0012)));
  }

  genes.forEach((gene) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "gene-button";
    button.dataset.gene = gene;
    button.textContent = gene;
    button.addEventListener("click", () => renderGene(gene));
    geneButtonsEl.appendChild(button);
  });

  canvas.addEventListener("pointerdown", onPointerDown);
  window.addEventListener("pointermove", onPointerMove);
  window.addEventListener("pointerup", onPointerUp);
  canvas.addEventListener("wheel", onWheel, { passive: false });
  window.addEventListener("resize", resize);

  resize();
  renderGene(genes[0]);
  animate();
})();
