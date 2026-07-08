import { PATCH_COLORS, faceColorIndex, vertexAngle, VERTEX_RATIO, N_PARTICLES } from "./constants.js";

// The simulation is genuinely 2D (particles live in a flat plane), so this is
// a pan/zoom canvas viewer rather than a 3D scene — rotating something flat
// out of its own plane would show nothing useful. Particle bodies are drawn
// exactly as codes/LatticeTools.py's `plot_2Dcolored_particle` does: six
// triangular wedges per hexagon, colored from a fixed 6-step viridis sample.

const PARTICLE_SIZE = 0.95; // matches particle_size used for the Triangular lattice in plot_a_system
const VERTEX_R = VERTEX_RATIO * PARTICLE_SIZE;

const ASYM_COLD = [0x39, 0x87, 0xe5]; // more stable than average
const ASYM_MID = [0x53, 0x51, 0x4c];
const ASYM_WARM = [0xe6, 0x67, 0x67]; // less stable than average

function lerpRgb(a, b, t) {
  return `rgb(${Math.round(a[0] + (b[0] - a[0]) * t)},${Math.round(a[1] + (b[1] - a[1]) * t)},${Math.round(a[2] + (b[2] - a[2]) * t)})`;
}

const VX = [0, 1, 2, 3, 4, 5].map((k) => Math.cos(vertexAngle(k)));
const VY = [0, 1, 2, 3, 4, 5].map((k) => Math.sin(vertexAngle(k)));

export function createViewer(container) {
  const canvas = document.createElement("canvas");
  container.appendChild(canvas);
  const ctx = canvas.getContext("2d");

  let cssW = 0, cssH = 0;
  let scale = 10; // px per world unit
  let center = { x: 0, y: 0 }; // world coords shown at canvas center
  let fitScale = 10;

  let positions = null; // Float32Array [x0,y0,x1,y1,...]
  let orientations = null;
  let asymmetry = null;
  let maxAbsAsym = 1;
  let mode = "patch"; // "patch" | "orientation" | "energy"
  let meta = null;

  function resize() {
    cssW = container.clientWidth;
    cssH = container.clientHeight;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = cssW * dpr;
    canvas.height = cssH * dpr;
    canvas.style.width = cssW + "px";
    canvas.style.height = cssH + "px";
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    draw();
  }
  new ResizeObserver(resize).observe(container);

  function worldToScreen(x, y) {
    return [cssW / 2 + (x - center.x) * scale, cssH / 2 - (y - center.y) * scale];
  }

  function particleUniformColor(i) {
    if (mode === "orientation") return PATCH_COLORS[orientations[i] - 1];
    if (mode === "energy") {
      const t = Math.max(-1, Math.min(1, asymmetry[i] / maxAbsAsym));
      return t < 0 ? lerpRgb(ASYM_MID, ASYM_COLD, -t) : lerpRgb(ASYM_MID, ASYM_WARM, t);
    }
    return null; // "patch" mode: per-wedge color, computed in draw()
  }

  function draw() {
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, cssW, cssH);
    if (!positions) return;

    for (let i = 0; i < N_PARTICLES; i++) {
      const x = positions[i * 2], y = positions[i * 2 + 1];
      const [sx, sy] = worldToScreen(x, y);
      if (sx < -30 || sx > cssW + 30 || sy < -30 || sy > cssH + 30) continue;

      const orientation = orientations[i];
      const uniform = particleUniformColor(i);

      for (let e = 0; e < 6; e++) {
        const x1 = sx + VX[e] * VERTEX_R * scale, y1 = sy - VY[e] * VERTEX_R * scale;
        const e2 = (e + 1) % 6;
        const x2 = sx + VX[e2] * VERTEX_R * scale, y2 = sy - VY[e2] * VERTEX_R * scale;
        ctx.beginPath();
        ctx.moveTo(sx, sy);
        ctx.lineTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.closePath();
        ctx.fillStyle = uniform || PATCH_COLORS[faceColorIndex(orientation, e)];
        ctx.fill();
      }
    }

    // periodic-boundary box frame -- particles wrap across this edge, so it's
    // not a decorative border, it's the actual simulation cell.
    const [bx0, by0] = worldToScreen(-meta.Lx / 2, -meta.Ly / 2);
    const [bx1, by1] = worldToScreen(meta.Lx / 2, meta.Ly / 2);
    ctx.strokeStyle = "#000000";
    ctx.lineWidth = 1;
    ctx.strokeRect(Math.min(bx0, bx1), Math.min(by0, by1), Math.abs(bx1 - bx0), Math.abs(by1 - by0));
  }

  function fit() {
    if (!meta) return;
    fitScale = Math.min(cssW / (meta.Lx * 1.12), cssH / (meta.Ly * 1.12));
    scale = fitScale;
    center = { x: 0, y: 0 };
  }

  function loadAssembly(record, entry) {
    meta = entry;
    positions = record.xy;
    orientations = record.orientation;
    asymmetry = record.asymmetry;
    maxAbsAsym = 1e-6;
    for (let i = 0; i < N_PARTICLES; i++) maxAbsAsym = Math.max(maxAbsAsym, Math.abs(asymmetry[i]));

    // recenter positions around the box center once, up front
    const cx = entry.Lx / 2, cy = entry.Ly / 2;
    for (let i = 0; i < N_PARTICLES; i++) {
      positions[i * 2] -= cx;
      positions[i * 2 + 1] -= cy;
    }

    resize();
    fit();
    draw();
  }

  function setMode(newMode) {
    mode = newMode;
    draw();
  }

  // --- pan / zoom interaction ---
  let dragging = false;
  let lastX = 0, lastY = 0;
  const activePointers = new Map();
  let pinchStartDist = 0, pinchStartScale = 1;

  canvas.style.touchAction = "none";
  canvas.style.cursor = "grab";

  canvas.addEventListener("pointerdown", (e) => {
    canvas.setPointerCapture(e.pointerId);
    activePointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
    if (activePointers.size === 1) {
      dragging = true;
      lastX = e.clientX;
      lastY = e.clientY;
      canvas.style.cursor = "grabbing";
    } else if (activePointers.size === 2) {
      dragging = false;
      const pts = [...activePointers.values()];
      pinchStartDist = Math.hypot(pts[0].x - pts[1].x, pts[0].y - pts[1].y) || 1;
      pinchStartScale = scale;
    }
  });

  canvas.addEventListener("pointermove", (e) => {
    if (!activePointers.has(e.pointerId)) return;
    activePointers.set(e.pointerId, { x: e.clientX, y: e.clientY });

    if (activePointers.size === 2) {
      const pts = [...activePointers.values()];
      const dist = Math.hypot(pts[0].x - pts[1].x, pts[0].y - pts[1].y) || 1;
      scale = Math.max(fitScale * 0.25, Math.min(fitScale * 25, pinchStartScale * (dist / pinchStartDist)));
      draw();
      return;
    }

    if (dragging) {
      const dx = e.clientX - lastX, dy = e.clientY - lastY;
      lastX = e.clientX; lastY = e.clientY;
      center.x -= dx / scale;
      center.y += dy / scale;
      draw();
    }
  });

  function endPointer(e) {
    activePointers.delete(e.pointerId);
    if (activePointers.size === 0) {
      dragging = false;
      canvas.style.cursor = "grab";
    }
  }
  canvas.addEventListener("pointerup", endPointer);
  canvas.addEventListener("pointercancel", endPointer);
  canvas.addEventListener("pointerleave", (e) => { if (activePointers.size <= 1) endPointer(e); });

  canvas.addEventListener("wheel", (e) => {
    e.preventDefault();
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left, my = e.clientY - rect.top;
    const worldX = center.x + (mx - cssW / 2) / scale;
    const worldY = center.y - (my - cssH / 2) / scale;
    const factor = Math.exp(-e.deltaY * 0.0015);
    scale = Math.max(fitScale * 0.25, Math.min(fitScale * 25, scale * factor));
    center.x = worldX - (mx - cssW / 2) / scale;
    center.y = worldY + (my - cssH / 2) / scale;
    draw();
  }, { passive: false });

  canvas.addEventListener("dblclick", () => { fit(); draw(); });

  return { loadAssembly, setMode };
}
