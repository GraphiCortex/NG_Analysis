
import React, {
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState
} from "react";

import ForceGraph3D from "react-force-graph-3d";
import * as THREE from "three";

const COLORS = {
  CORE: "#f0a579",
  INTERFACE: "#5bc7e7",
  BRIDGE: "#fc66c0",
  EXTERNAL: "#a18cff",
  GOLD: "#ffcc65",
  ROUTE: "#92f1be"
};

const STARTERS = [
  "RAD51", "TIMELESS", "TOPBP1",
  "CLSPN", "RPA2", "FAS"
];


const GENE_ALIASES = {
  H2AX: "H2AFX", NOXA: "PMAIP1", PUMA: "BBC3",
  DAI: "ZBP1", DLM1: "ZBP1", HGFR: "MET"
};

const EMPTY_PATH = [];

function matchingGenes(value, names) {
  const query = value.trim().toUpperCase();
  if (!query) return [];

  const matches = new Set();

  if (
    GENE_ALIASES[query] &&
    names.includes(GENE_ALIASES[query])
  ) {
    matches.add(GENE_ALIASES[query]);
  }

  for (const [alias, gene] of Object.entries(GENE_ALIASES)) {
    if (alias.startsWith(query) && names.includes(gene))
      matches.add(gene);
  }

  for (const name of names) {
    if (name.includes(query))
      matches.add(name);
  }

  return [...matches].sort((a, b) => {
    const rank = name =>
      name === query ? 0 :
      name.startsWith(query) ? 1 : 2;

    return rank(a) - rank(b) || a.localeCompare(b);
  }).slice(0, 8);
}

const PAIR_REVIEW = {
  "RAD51|CASP3":
    "Reviewed cleavage evidence points from CASP3 to RAD51; the database direction is questionable.",

  "RAD51|BRCA2":
    "BRCA2–RAD51 is an established partnership; RAD51 activating BRCA2 is not established.",

  "BRCA2|SMAD3":
    "A physical/functional interaction is reported, not necessarily a directed signaling cascade.",

  "TIMELESS|ATM":
    "Evidence supports a checkpoint role, but not direct activation of ATM by TIMELESS.",

  "ATM|FAS":
    "Evidence is context-specific and has not been verified in the Fong neural stem-cell system.",

  "ATM|TNFRSF10B":
    "The direction or sign of regulation is context-dependent in the reviewed literature.",

  "MYC|FAS":
    "The reviewed transcriptional evidence conflicts with the encoded inhibitory sign.",

  "MYC|TNFRSF10B":
    "Supporting effects were reported in cancer models, not verified in the Fong cells.",

  "TP53|PYCARD":
    "p53-dependent expression has been reported; direct promoter regulation is a separate question.",

  "PTK2B|PYCARD":
    "ASC phosphorylation is supported in inflammasome studies; neural relevance is not established.",

  "TOPBP1|TP53":
    "The exported database edge has no original reference tokens.",

  "TNF|MLKL":
    "This arrow summarizes a multistep necroptotic process rather than a direct biochemical action."
};

const edgeKey = (a, b) => `${a}|${b}`;

const signText = x => ({
  stimulation: "Stimulatory (database)",
  inhibition: "Inhibitory (database)",
  unsigned: "Effect unknown",
  ambiguous: "Effect uncertain"
}[x] || "Effect not specified");

const getColor = n =>
  n.tier?.includes("EXTERNAL")
    ? COLORS.EXTERNAL
    : n.f4_cross_module_bridge
      ? COLORS.BRIDGE
      : (COLORS[n.tier] || COLORS.EXTERNAL);

/*
 * Original visual quality, optimized internally.
 *
 * Shared geometry and materials reduce GPU allocations.
 * High-resolution labels are generated only once per gene.
 */

const SPHERE = new THREE.SphereGeometry(1, 22, 16);
const OCTAHEDRON = new THREE.OctahedronGeometry(1, 0);
const HALO_SPHERE = new THREE.SphereGeometry(1, 14, 10);
const RING = new THREE.TorusGeometry(
  1.2, 0.075, 7, 35
);

const BODY_MATERIALS = new Map();
const HALO_MATERIALS = new Map();

const RING_MATERIAL = new THREE.MeshBasicMaterial({
  color: COLORS.GOLD
});

const PLATE_MATERIALS = new Map();

function bodyMaterial(color, highlighted = false) {
  const key = `${color}:${highlighted}`;

  if (!BODY_MATERIALS.has(key)) {
    BODY_MATERIALS.set(
      key,
      new THREE.MeshStandardMaterial({
        color,
        emissive: color,
        emissiveIntensity: highlighted ? 0.8 : 0.32,
        metalness: 0.14,
        roughness: 0.46
      })
    );
  }

  return BODY_MATERIALS.get(key);
}

function haloMaterial(color) {
  if (!HALO_MATERIALS.has(color)) {
    HALO_MATERIALS.set(
      color,
      new THREE.MeshBasicMaterial({
        color,
        transparent: true,
        opacity: 0.10,
        depthWrite: false
      })
    );
  }

  return HALO_MATERIALS.get(color);
}

function nameplate(label, color) {
  const key = `${label}|${color}`;

  let material = PLATE_MATERIALS.get(key);

  if (!material) {
    const canvas = document.createElement("canvas");

    canvas.width = 512;
    canvas.height = 144;

    const ctx = canvas.getContext("2d");

    ctx.fillStyle = "rgba(9,27,43,0.90)";
    ctx.strokeStyle = color;
    ctx.lineWidth = 5;

    ctx.beginPath();
    ctx.roundRect(7, 9, 498, 126, 45);
    ctx.fill();
    ctx.stroke();

    ctx.fillStyle = "#f7fbff";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.font = "700 63px system-ui, sans-serif";

    ctx.fillText(label, 256, 72, 456);

    const texture = new THREE.CanvasTexture(canvas);
    texture.colorSpace = THREE.SRGBColorSpace;

    material = new THREE.SpriteMaterial({
      map: texture,
      transparent: true,
      depthTest: false,
      depthWrite: false
    });

    PLATE_MATERIALS.set(key, material);
  }

  const sprite = new THREE.Sprite(material);

  sprite.scale.set(12.4, 3.48, 1);
  sprite.position.set(0, -6.8, 0);
  sprite.renderOrder = 12;

  return sprite;
}

function Note({ children }) {
  return <p className="micro-note">{children}</p>;
}


function GeneAutocomplete({
  value, onChange, onSelect, names, placeholder
}) {
  const [open, setOpen] = useState(false);
  const [active, setActive] = useState(0);

  const results = useMemo(
    () => matchingGenes(value, names),
    [value, names]
  );

  const select = gene => {
    onSelect(gene);
    setOpen(false);
    setActive(0);
  };

  return (
    <div className="gene-search" style={{ minWidth: 0 }}>
      <input
        className="field"
        type="search"
        autoComplete="off"
        spellCheck={false}
        value={value}
        placeholder={placeholder}
        aria-label="Gene symbol"
        onFocus={() => setOpen(true)}
        onBlur={() => setOpen(false)}
        onChange={e => {
          onChange(e.target.value.toUpperCase());
          setOpen(true);
          setActive(0);
        }}
        onKeyDown={e => {
          if (e.key === "Escape")
            setOpen(false);

          if (e.key === "ArrowDown" && results.length) {
            e.preventDefault();
            setActive(i => (i + 1) % results.length);
          }

          if (e.key === "ArrowUp" && results.length) {
            e.preventDefault();
            setActive(i =>
              (i - 1 + results.length) % results.length
            );
          }

          if (e.key === "Enter" && results.length) {
            e.preventDefault();
            select(
              results[Math.min(active, results.length - 1)]
            );
          }
        }}
      />

      {open && value.trim() && (
        <div
          className="suggestions"
          role="listbox"
          style={{ maxHeight: 220, overflowY: "auto" }}
        >
          {results.length ? results.map((gene, i) => (
            <button
              type="button"
              key={gene}
              role="option"
              aria-selected={i === active}
              style={
                i === active
                  ? { background: "#295775" }
                  : undefined
              }
              onMouseDown={e => e.preventDefault()}
              onClick={() => select(gene)}
            >
              {gene}
              {GENE_ALIASES[value.trim().toUpperCase()] ===
                gene &&
               gene !== value.trim().toUpperCase()
                ? " (alias)" : ""}
            </button>
          )) : (
            <div style={{
              padding: 11,
              background: "#173149",
              color: "#b7cddc",
              fontSize: 11
            }}>
              No matching gene in the exported atlas.
            </div>
          )}
        </div>
      )}
    </div>
  );
}

// This is deliberately a preview search over the
// 393 saved edges, not the full OmniPath network.
function searchSavedPaths(atlas, source, target, key) {
  const bounds = atlas.metadata.settings[key];

  const all = new Map(
    atlas.nodes.map(node => [node.id, node])
  );

  if (
    !bounds ||
    !all.has(source) ||
    !all.has(target) ||
    source === target
  ) {
    return [];
  }

  const next = new Map();

  for (const edge of atlas.directed_edges) {
    if (!next.has(edge.source)) {
      next.set(edge.source, []);
    }

    next.get(edge.source).push(edge);
  }

  const queue = [{
    nodes: [source],
    edges: [],
    outside: 0
  }];

  const found = [];
  let visited = 0;

  for (
    let head = 0;
    head < queue.length &&
    visited < 20000 &&
    found.length < 12;
    head++
  ) {
    visited++;

    const route = queue[head];
    const last = route.nodes.at(-1);

    if (last === target) {
      found.push(route);
      continue;
    }

    if (route.edges.length >= bounds.max_edges) {
      continue;
    }

    for (const edge of next.get(last) || []) {
      if (route.nodes.includes(edge.target)) {
        continue;
      }

      const extra = route.outside + Number(
        edge.target !== target &&
        all.get(edge.target)?.tier === "F5_EXTERNAL"
      );

      if (extra > bounds.max_external) {
        continue;
      }

      queue.push({
        nodes: [...route.nodes, edge.target],
        edges: [...route.edges, edge],
        outside: extra
      });
    }
  }

  return found;
}

export default function App() {
  const [atlas, setAtlas] = useState(null);
  const [research, setResearch] = useState(null);
  const [error, setError] = useState("");

  const [mode, setMode] = useState("focus");
  const [focus, setFocus] = useState("RAD51");
  const [selected, setSelected] = useState("RAD51");
  const [inspector, setInspector] = useState(true);

  const [geneQuery, setGeneQuery] = useState("");

  // Neighborhood customization.
  const [scope, setScope] = useState("all");
  const [picked, setPicked] = useState([]);
  const [neighborQuery, setNeighborQuery] = useState("");
  const [showBetween, setShowBetween] = useState(false);

  // Network layers.
  const [showDirect, setShowDirect] = useState(false);
  const [showSaved, setShowSaved] = useState(false);
  const [showExternal, setShowExternal] = useState(false);

  const [moduleId, setModuleId] = useState("ALL");
  const [bridgesOnly, setBridgesOnly] = useState(false);
  const [networkFocus, setNetworkFocus] = useState(null);

  // Path explorer.
  const [pathKind, setPathKind] = useState("examples");
  const [source, setSource] = useState("RAD51");
  const [target, setTarget] = useState("FAS");
  const [from, setFrom] = useState("E2F1");
  const [to, setTo] = useState("FAS");
  const [setting, setSetting] = useState("E4_X1");
  const [routeIndex, setRouteIndex] = useState(0);

  const [dimensions, setDimensions] = useState({
    width: 1100,
    height: 700
  });

  const canvasRef = useRef(null);
  const graphRef = useRef(null);

  const pendingAtlasClick = useRef({
    gene: null,
    time: 0,
    timer: null
  });

  useEffect(() => () => {
    if (pendingAtlasClick.current.timer) {
      clearTimeout(pendingAtlasClick.current.timer);
    }
  }, []);

  // =====================================================
  // LOAD DATA
  // =====================================================

  useEffect(() => {
    const controller = new AbortController();

    const read = filename =>
      fetch(
        `${import.meta.env.BASE_URL}data/${filename}`,
        { signal: controller.signal }
      ).then(response => {
        if (!response.ok) {
          throw Error(
            `${filename}: HTTP ${response.status}`
          );
        }

        return response.json();
      });

    Promise.all([
      read("F8D_viewer_data.json"),
      read("F9_overlay.json")
    ])
      .then(([base, extra]) => {
        if (
          base.nodes?.length !== 343 ||
          base.string_edges?.length !== 1139 ||
          base.directed_edges?.length !== 393 ||
          extra.directEdges?.length !== 963 ||
          extra.bypassRoutes?.length !== 25
        ) {
          throw Error(
            "The network exports do not match this website version."
          );
        }

        setAtlas(base);
        setResearch(extra);
      })
      .catch(e => {
        if (e.name !== "AbortError") {
          setError(e.message);
        }
      });

    return () => controller.abort();
  }, []);

  useEffect(() => {
    if (!canvasRef.current) return;

    const observer = new ResizeObserver(([entry]) => {
      if (entry) {
        setDimensions({
          width: entry.contentRect.width,
          height: entry.contentRect.height
        });
      }
    });

    observer.observe(canvasRef.current);

    return () => observer.disconnect();
  }, [atlas]);

  // =====================================================
  // INDEXES
  // =====================================================

  const byId = useMemo(
    () => new Map(
      (atlas?.nodes || []).map(g => [g.id, g])
    ),
    [atlas]
  );

  const modules = useMemo(
    () => new Map(
      (atlas?.modules || []).map(m => [
        m.id,
        m.label
      ])
    ),
    [atlas]
  );

  const geneNames = useMemo(
    () => [...byId.keys()].sort(),
    [byId]
  );

  // Combine the three recorded edge sources.
  const allEdges = useMemo(() => {
    if (!atlas || !research) return [];

    return [
      ...atlas.string_edges.map(e => ({
        source: e.source,
        target: e.target,
        kind: "association",
        score: e.score_1000
      })),

      ...research.directEdges.map(e => ({
        source: e.source,
        target: e.target,
        kind: "direct",
        effect: e.effect
      })),

      ...atlas.directed_edges.map(e => ({
        source: e.source,
        target: e.target,
        kind: "saved",
        effect: e.effect_label
      }))
    ];
  }, [atlas, research]);

  // =====================================================
  // DIRECT NEIGHBOR OPTIONS
  // =====================================================

  // The picker uses all available edge types,
  // irrespective of the main network-layer toggles.
  const neighbors = useMemo(() => {
    const options = new Map();

    for (const edge of allEdges) {
      if (
        edge.source !== focus &&
        edge.target !== focus
      ) {
        continue;
      }

      const other =
        edge.source === focus
          ? edge.target
          : edge.source;

      if (!options.has(other)) {
        options.set(other, {
          gene: other,
          types: new Set(),
          score: 0
        });
      }

      const option = options.get(other);

      option.types.add(edge.kind);

      option.score = Math.max(
        option.score,
        edge.score || 0
      );
    }

    return [...options.values()].sort((a, b) =>
      Number(b.types.has("direct")) -
      Number(a.types.has("direct")) ||
      b.score - a.score ||
      a.gene.localeCompare(b.gene)
    );
  }, [allEdges, focus]);

  const optionsShown = useMemo(
    () => neighbors.filter(
      n => n.gene.includes(
        neighborQuery.trim().toUpperCase()
      )
    ).slice(0, 60),
    [neighbors, neighborQuery]
  );

  const chosenSet = useMemo(
    () => new Set(picked),
    [picked]
  );


  const pickFocus = useCallback(gene => {
    setFocus(gene);
    setSelected(gene);
    setMode("focus");
    setInspector(true);
    setScope("all");
    setPicked([]);
    setNeighborQuery("");
    setGeneQuery("");

    const node = byId.get(gene);
    const hasString = atlas?.string_edges.some(
      e => e.source === gene || e.target === gene
    ) ?? false;

    if (
      node?.tier?.includes("EXTERNAL") ||
      !hasString
    ) {
      setShowExternal(true);
      setShowDirect(true);
      setShowSaved(true);
    }
  }, [byId, atlas]);

  // Click once to inspect a gene.
  // Double-click to display its immediate neighborhood.
  // Single click: inspect.
  // Double click: local view INSIDE the selected network.
  const handleNodeClick = useCallback(node => {
    const gene = node.id;

    if (mode !== "atlas") {
      if (mode === "focus") {
        pickFocus(gene);
      } else {
        setSelected(gene);
        setInspector(true);
      }
      return;
    }

    const last = pendingAtlasClick.current;
    const now = Date.now();

    if (
      last.gene === gene &&
      now - last.time < 420
    ) {
      if (last.timer) clearTimeout(last.timer);

      pendingAtlasClick.current = {
        gene: null, time: 0, timer: null
      };

      // Stay in Network mode, retaining the process.
      setNetworkFocus(gene);
      setSelected(gene);
      setInspector(true);
      return;
    }

    if (last.timer) clearTimeout(last.timer);

    const timer = setTimeout(() => {
      setSelected(gene);
      setInspector(true);

      if (pendingAtlasClick.current.timer === timer) {
        pendingAtlasClick.current = {
          gene: null, time: 0, timer: null
        };
      }
    }, 420);

    pendingAtlasClick.current = {
      gene, time: now, timer
    };
  }, [mode, pickFocus]);

  const togglePick = useCallback(
    gene => setPicked(current =>
      current.includes(gene)
        ? current.filter(g => g !== gene)
        : [...current, gene]
    ),
    []
  );

  // =====================================================
  // ROUTE DATA
  // =====================================================

  const availableRouteSources = useMemo(
    () => [...new Set(
      (research?.bypassRoutes || [])
        .map(r => r.source)
    )],
    [research]
  );

  const availableRouteTargets = useMemo(
    () => [...new Set(
      (research?.bypassRoutes || [])
        .filter(r => r.source === source)
        .map(r => r.target)
    )],
    [research, source]
  );

  const examples = useMemo(
    () => (research?.bypassRoutes || [])
      .filter(r =>
        r.source === source &&
        r.target === target
      ),
    [research, source, target]
  );

  const legacyPaths = useMemo(
    () => atlas
      ? searchSavedPaths(
          atlas,
          from.trim().toUpperCase(),
          to.trim().toUpperCase(),
          setting
        )
      : [],
    [atlas, from, to, setting]
  );

  const paths =
    pathKind === "examples"
      ? examples
      : legacyPaths;

  const activePath = paths[
    Math.min(routeIndex, paths.length - 1)
  ];

  const pathNodes = activePath?.nodes || EMPTY_PATH;

  const bridgeSet = useMemo(
    () => new Set(
      (atlas?.nodes || [])
        .filter(n =>
          n.f4_cross_module_bridge &&
          (
            moduleId === "ALL" ||
            n.modules.includes(moduleId)
          )
        )
        .map(n => n.id)
    ),
    [atlas, moduleId]
  );

  // =====================================================
  // VISUAL EDGE SELECTION
  // =====================================================

  // Do not render overlapping network layers
  // multiple times for the same connection.
  const visualEdges = useMemo(() => {
    const chosen = new Map();

    for (const edge of allEdges) {
      if (
        scope === "choose" &&
        mode === "focus"
      ) {
        const focal =
          edge.source === focus ||
          edge.target === focus;

        if (
          focal &&
          !chosenSet.has(
            edge.source === focus
              ? edge.target
              : edge.source
          )
        ) {
          continue;
        }

        if (!focal && !showBetween) continue;

        if (
          !focal &&
          (
            !chosenSet.has(edge.source) ||
            !chosenSet.has(edge.target)
          )
        ) {
          continue;
        }
      } else if (
        edge.kind === "direct" &&
        !showDirect
      ) {
        continue;
      } else if (
        edge.kind === "saved" &&
        !showSaved
      ) {
        continue;
      }

      const isUndirected =
        edge.kind === "association";

      const key = isUndirected
        ? [edge.source, edge.target]
            .sort()
            .join("|")
        : edgeKey(
            edge.source,
            edge.target
          );

      const previous = chosen.get(key);

      const rank = {
        association: 1,
        saved: 2,
        direct: 3
      };

      if (
        !previous ||
        rank[edge.kind] >
        rank[previous.kind]
      ) {
        chosen.set(key, edge);
      }
    }

    return [...chosen.values()];
  }, [
    allEdges,
    mode,
    scope,
    focus,
    chosenSet,
    showBetween,
    showDirect,
    showSaved
  ]);

  // =====================================================
  // 3D GRAPH DATA
  // =====================================================

  const graphData = useMemo(() => {
    if (!atlas || !research) {
      return {
        nodes: [],
        links: []
      };
    }

    // -------------------------------------------
    // PATH VIEW
    // -------------------------------------------

    if (mode === "route") {
      if (!activePath) {
        return {
          nodes: [],
          links: []
        };
      }

      const nodes = pathNodes.map(id => ({
        ...(
          byId.get(id) || {
            id,
            tier: "F9_EXTERNAL",
            modules: []
          }
        ),
        showTag: true
      }));

      const links = pathNodes
        .slice(1)
        .map((id, i) => {
          const a = pathNodes[i];

          const evidence =
            research.edgeEvidence[
              edgeKey(a, id)
            ];

          const original =
            pathKind === "saved"
              ? activePath.edges[i]
              : null;

          return {
            id: `path${i}`,
            source: a,
            target: id,
            kind: "route",
            effect:
              evidence?.effect ||
              original?.effect_label ||
              "unsigned"
          };
        });

      return { nodes, links };
    }

    // -------------------------------------------
    // GENE NEIGHBORHOOD
    // -------------------------------------------

    if (mode === "focus") {
      const focalEdges = visualEdges.filter(
        e =>
          e.source === focus ||
          e.target === focus
      );

      const adjacent = new Set([focus]);

      for (const edge of focalEdges) {
        adjacent.add(edge.source);
        adjacent.add(edge.target);
      }

      const ids =
        scope === "choose"
          ? new Set([focus, ...picked])
          : adjacent;

      const links = visualEdges.filter(e => {
        if (
          !ids.has(e.source) ||
          !ids.has(e.target)
        ) {
          return false;
        }

        return (
          (
            showBetween &&
            scope === "choose"
          ) ||
          e.source === focus ||
          e.target === focus
        );
      });

      const nodes = [...ids]
        .map(id =>
          byId.get(id) || {
            id,
            tier: "F9_EXTERNAL",
            modules: []
          }
        )
        .filter(n =>
          showExternal ||
          !n.tier?.includes("EXTERNAL") ||
          scope === "choose"
        )
        .map(n => ({
          ...n,
          showTag:
            ids.size <= 16 ||
            n.id === focus
        }));

      const retained = new Set(
        nodes.map(n => n.id)
      );

      return {
        nodes,
        links: links
          .filter(e =>
            retained.has(e.source) &&
            retained.has(e.target)
          )
          .map((e, i) => ({
            ...e,
            id:
              `${e.kind}:${e.source}:${e.target}:${i}`
          }))
      };
    }

    // -------------------------------------------
    // WHOLE NETWORK (filtered local drilldown)
    // -------------------------------------------
    // Network view: apply biological filtering FIRST.
    let networkNodes;
    let networkLinks;

    if (bridgesOnly) {
      const connectorLinks = atlas.string_edges.filter(
        edge =>
          (
            bridgeSet.has(edge.source) &&
            byId.get(edge.target)?.tier === "CORE"
          ) ||
          (
            bridgeSet.has(edge.target) &&
            byId.get(edge.source)?.tier === "CORE"
          )
      );

      const allowed = new Set(bridgeSet);

      for (const edge of connectorLinks) {
        allowed.add(edge.source);
        allowed.add(edge.target);
      }

      networkNodes = atlas.nodes.filter(
        n => allowed.has(n.id)
      );

      networkLinks = connectorLinks.map(
        (edge, i) => ({
          ...edge,
          id: `connector:${i}`,
          kind: "association"
        })
      );
    } else {
      // Only genes annotated with the selected
      // biological process are included.
      networkNodes = atlas.nodes.filter(
        n =>
          (
            showExternal ||
            n.tier !== "F5_EXTERNAL"
          ) &&
          (
            moduleId === "ALL" ||
            n.modules.includes(moduleId)
          )
      );

      const allowed = new Set(
        networkNodes.map(n => n.id)
      );

      networkLinks = visualEdges.filter(
        edge =>
          allowed.has(edge.source) &&
          allowed.has(edge.target)
      ).map((edge, i) => ({
        ...edge,
        id: `atlas:${edge.kind}:${i}`
      }));
    }

    // Double-click: ONE-HOP neighborhood within
    // the already-filtered graph.
    if (
      networkFocus &&
      networkNodes.some(n => n.id === networkFocus)
    ) {
      const local = new Set([networkFocus]);

      for (const edge of networkLinks) {
        if (edge.source === networkFocus) {
          local.add(edge.target);
        }

        if (edge.target === networkFocus) {
          local.add(edge.source);
        }
      }

      networkNodes = networkNodes.filter(
        n => local.has(n.id)
      );

      // Retain existing interactions among the
      // local neighbors, provided those interactions
      // belong to the selected network.
      networkLinks = networkLinks.filter(
        edge =>
          local.has(edge.source) &&
          local.has(edge.target)
      );
    }

    return {
      nodes: networkNodes.map(n => ({
        ...n,
        showTag:
          !!networkFocus &&
          (
            networkNodes.length <= 16 ||
            n.id === networkFocus
          )
      })),
      links: networkLinks
    };

  }, [
    atlas,
    research,
    mode,
    focus,
    scope,
    picked,
    visualEdges,
    showExternal,
    showBetween,
    moduleId,
    bridgesOnly,
    networkFocus,
    bridgeSet,
    byId,
    activePath,
    pathNodes,
    pathKind
  ]);

  // =====================================================
  // LIGHTWEIGHT THREE.JS NODES
  // =====================================================

  
const makeNode = useCallback(n => {
  const group = new THREE.Group();
  const color = getColor(n);

  const external = n.tier?.includes("EXTERNAL");

  const radius = n.f4_cross_module_bridge
    ? 5.7
    : external
      ? 4.8
      : 4.5;

  // Restore detailed spherical geometry and the
  // emissive material used by the original website.
  const mesh = new THREE.Mesh(
    external ? OCTAHEDRON : SPHERE,
    bodyMaterial(color, n.id === selected)
  );

  mesh.scale.setScalar(radius);
  group.add(mesh);

  // Restore the translucent magenta bridge halo.
  if (n.f4_cross_module_bridge) {
    const halo = new THREE.Mesh(
      HALO_SPHERE,
      haloMaterial(color)
    );

    halo.scale.setScalar(radius * 1.16);

    group.add(halo);
  }

  // Restore the gold overlap ring.
  if (n.shared37) {
    const ring = new THREE.Mesh(
      RING,
      RING_MATERIAL
    );

    ring.scale.setScalar(radius);
    group.add(ring);
  }

  // Restore the rounded, high-resolution nameplates.
  // Avoid drawing unnecessary labels in large views.
  if (
    n.showTag ||
    n.id === selected ||
    (mode === "focus" && n.id === focus)
  ) {
    group.add(nameplate(n.id, color));
  }

  return group;
}, [selected, focus, mode]);


  // =====================================================
  // CAMERA
  // =====================================================

  const fit = useCallback(() => {
    const g = graphRef.current;

    if (!g || !graphData.nodes.length) {
      return;
    }

    if (graphData.nodes.length === 1) {
      g.cameraPosition(
        { x: 0, y: 0, z: 100 },
        { x: 0, y: 0, z: 0 },
        400
      );
    } else {
      g.zoomToFit(450, 95);
    }
  }, [graphData.nodes.length]);

  useEffect(() => {
    if (!atlas || !research) {
      return;
    }

    const t = setTimeout(fit, 350);

    return () => clearTimeout(t);
  }, [
    atlas,
    research,
    graphData,
    fit
  ]);

  // =====================================================
  // INSPECTOR DATA
  // =====================================================

  const details =
    byId.get(selected) ||
    (
      selected
        ? {
            id: selected,
            tier: "F9_EXTERNAL",
            modules: []
          }
        : null
    );

  const expression =
    research?.genes?.[selected];

  const allSelectedEdges = useMemo(
    () => allEdges.filter(e =>
      e.source === selected ||
      e.target === selected
    ),
    [allEdges, selected]
  );

  const getRouteDetails =
    activePath &&
    mode === "route" &&
    pathKind === "examples"
      ? pathNodes.slice(1).map((id, i) => ({
          source: pathNodes[i],
          target: id,

          evidence:
            research?.edgeEvidence[
              edgeKey(pathNodes[i], id)
            ],

          warning:
            PAIR_REVIEW[
              edgeKey(pathNodes[i], id)
            ]
        }))
      : [];

  if (error) {
    return (
      <div className="loading">
        Unable to open network: {error}
      </div>
    );
  }

  if (!atlas || !research) {
    return (
      <div className="loading">
        Opening the gene atlas…
      </div>
    );
  }

  // =====================================================
  // RENDER
  // =====================================================

  return (
    <div className="shell">
      <header className="top">
        <div className="brand">
          <span className="brand-mark">
            R·N
          </span>

          <div>
            <div className="eyebrow">
              NOËL LAB · RESEARCH ATLAS
            </div>

            <h1>
              Neural cell fate{" "}
              <em>· 3D gene network</em>
            </h1>
          </div>
        </div>

        <div className="top-right">
          Mouse neural stem cells · Fong et al. (2022)
        </div>
      </header>

      <main
        className="stage"
        ref={canvasRef}
      >
        <ForceGraph3D
          ref={graphRef}
          width={dimensions.width}
          height={dimensions.height}
          graphData={graphData}
          backgroundColor="#07111d"
          showNavInfo={false}

          nodeThreeObject={makeNode}

          nodeLabel={n =>
            `${n.id} · ${
              n.tier?.includes("EXTERNAL")
                ? "Additional protein"
                : n.tier === "CORE"
                  ? "Cell-fate gene"
                  : "Related-process gene"
            }`
          }

          onNodeClick={handleNodeClick}

          onBackgroundClick={() =>
            setInspector(false)
          }

          linkColor={e =>
            e.kind === "route"
              ? COLORS.ROUTE
              : e.kind === "direct"
                ? "#87bae8"
                : e.kind === "saved"
                  ? "#a18cff"
                  : "#496479"
          }

          linkOpacity={0.47}


          linkWidth={e =>
            graphData.links.length > 160
              ? 0
              : e.kind === "route"
                ? 1.8
                : e.kind === "association"
                  ? 0.42
                  : 0.9
          }

          // SCOPED_NETWORK_ARROWS
          // Show direction in small biological networks,
          // but avoid excessive arrow geometry in the atlas.
          linkDirectionalArrowLength={e => {
            if (e.kind === "association") return 0;

            if (mode === "route") return 2.8;

            if (
              mode === "focus" &&
              graphData.links.length <= 60
            ) return 2.4;

            if (
              mode === "atlas" &&
              (moduleId !== "ALL" || networkFocus !== null) &&
              graphData.links.length <= 180
            ) return 5.0;

            return 0;
          }}

          // DIRECTIONAL_VISIBILITY_FIX
          linkDirectionalArrowRelPos={0.79}
          linkDirectionalArrowResolution={8}

          linkDirectionalArrowColor={e =>
            e.effect === "inhibition" ? "#ff9bac" :
            e.kind === "direct" ? "#a8f1ff" :
            e.kind === "saved" ? "#d4b3ff" :
            COLORS.ROUTE
          }

          // Direction is emphasized only in small
          // process-specific networks.
          linkDirectionalParticles={e =>
            mode === "atlas" &&
            (moduleId !== "ALL" || networkFocus !== null) &&
            graphData.links.length <= 65 &&
            (e.kind === "direct" || e.kind === "saved")
              ? 1 : 0
          }

          linkDirectionalParticleWidth={1.6}
          linkDirectionalParticleSpeed={0.006}
          linkDirectionalParticleResolution={6}

          linkDirectionalParticleColor={e =>
            e.effect === "inhibition" ? "#ff9bac" :
            e.kind === "saved" ? "#d4b3ff" :
            "#a8f1ff"
          }

          // Separate overlapping directed records.
          linkCurvature={e =>
            mode === "route" ? 0.06 :
            mode === "atlas" &&
            e.kind !== "association" &&
            graphData.links.length <= 80 &&
            (moduleId !== "ALL" || networkFocus !== null)
              ? 0.15 : 0
          }

          // Hover to inspect the connection type.
          linkLabel={e => {
            const a = typeof e.source === "object"
              ? e.source.id : e.source;
            const b = typeof e.target === "object"
              ? e.target.id : e.target;

            if (e.kind === "association") {
              return `${a} -- ${b} | STRING association`;
            }

            const kind =
              e.kind === "direct"
                ? "Directed database interaction"
                : e.kind === "saved"
                  ? "Earlier illustrative interaction"
                  : "Selected route";

            return `${a} -> ${b} | ${kind} | ${signText(e.effect)}`;
          }}

          nodeResolution={10}
          cooldownTicks={65}
          d3AlphaDecay={0.045}

          enableNodeDrag={
            graphData.nodes.length < 120
          }
        />

        {/* LEFT PANEL */}

        <aside className="control-panel">
          <div className="panel-kicker">
            EXPLORE
          </div>

          <details className="compact-info">
            <summary>About the map</summary>
            <p>
              Explore gene-expression changes in mouse
              neural stem cells alongside published
              interaction networks.
            </p>
          </details>

          <div className="mode-switch">
            {[
              ["focus", "Gene"],
              ["atlas", "Network"],
              ["route", "Routes"]
            ].map(([id, text]) => (
              <button
                key={id}
                className={
                  mode === id ? "current" : ""
                }
                onClick={() => {
                  setMode(id);
                  setInspector(false);
                  setRouteIndex(0);
                  setNetworkFocus(null);
                }}
              >
                {text}
              </button>
            ))}
          </div>

          {/* GENE EXPLORER */}

          {mode === "focus" && (
            <>
              <label className="field-label">
                Find a gene
              </label>


              <GeneAutocomplete
                value={geneQuery}
                onChange={setGeneQuery}
                onSelect={pickFocus}
                names={geneNames}
                placeholder="e.g. RAD51, CHEK1, FAS"
              />

              <div className="quick-genes">
                {STARTERS.map(gene => (
                  <button
                    key={gene}
                    className={
                      focus === gene
                        ? "chosen"
                        : ""
                    }
                    onClick={() =>
                      pickFocus(gene)
                    }
                  >
                    {gene}
                  </button>
                ))}
              </div>

              {/* CUSTOM SUBNEIGHBORHOOD */}

              <div className="subset-panel">
                <strong>
                  Explore {focus}'s connections
                </strong>

                <div className="subset-switch">
                  <button
                    className={
                      scope === "all"
                        ? "active"
                        : ""
                    }
                    onClick={() =>
                      setScope("all")
                    }
                  >
                    All neighbors
                  </button>

                  <button
                    className={
                      scope === "choose"
                        ? "active"
                        : ""
                    }
                    onClick={() =>
                      setScope("choose")
                    }
                  >
                    Choose genes
                  </button>
                </div>

                {scope === "choose" ? (
                  <>
                    <Note>
                      Choose exactly which
                      genes to compare with{" "}
                      {focus}. Only those
                      genes and their recorded
                      links will appear.
                    </Note>

                    <input
                      className="field"
                      value={neighborQuery}
                      placeholder="Search neighbors…"
                      onChange={e =>
                        setNeighborQuery(
                          e.target.value
                        )
                      }
                    />

                    <div className="subset-actions">
                      <button
                        onClick={() =>
                          setPicked(
                            optionsShown
                              .slice(0, 8)
                              .map(n => n.gene)
                          )
                        }
                      >
                        Pick first 8
                      </button>

                      <button
                        onClick={() =>
                          setPicked([])
                        }
                      >
                        Clear
                      </button>

                      <span>
                        {picked.length} selected
                      </span>
                    </div>

                    <div className="subset-options">
                      {optionsShown.map(item => (
                        <label
                          key={item.gene}
                          className="subset-option"
                        >
                          <input
                            type="checkbox"
                            checked={
                              chosenSet.has(
                                item.gene
                              )
                            }
                            onChange={() =>
                              togglePick(item.gene)
                            }
                          />

                          <span>
                            <b>{item.gene}</b>

                            <small>
                              {item.types.has("direct")
                                ? "Directed"
                                : item.types.has("saved")
                                  ? "Earlier route"
                                  : "STRING association"}
                            </small>
                          </span>
                        </label>
                      ))}

                      {!optionsShown.length && (
                        <Note>
                          No matching direct
                          neighbors in the
                          exported data.
                        </Note>
                      )}
                    </div>

                    <label className="mini-check">
                      <input
                        type="checkbox"
                        checked={showBetween}
                        onChange={e =>
                          setShowBetween(
                            e.target.checked
                          )
                        }
                      />

                      Also show connections
                      among chosen genes
                    </label>

                    {!picked.length && (
                      <Note>
                        Select a gene above
                        to build your subnetwork.
                      </Note>
                    )}
                  </>
                ) : (
                  <Note>
                    {neighbors.length} recorded
                    neighbors across all layers.
                    Only links touching {focus}
                    are drawn to reduce clutter.
                  </Note>
                )}
              </div>

              {scope === "all" && (
                <details className="layers">
                  <summary>
                    Connection layers
                  </summary>

                  <label>
                    <input
                      type="checkbox"
                      checked={showDirect}
                      onChange={e =>
                        setShowDirect(
                          e.target.checked
                        )
                      }
                    />
                    Directed database interactions
                  </label>

                  <label>
                    <input
                      type="checkbox"
                      checked={showSaved}
                      onChange={e =>
                        setShowSaved(
                          e.target.checked
                        )
                      }
                    />
                    Earlier illustrative interactions
                  </label>

                  <label>
                    <input
                      type="checkbox"
                      checked={showExternal}
                      onChange={e =>
                        setShowExternal(
                          e.target.checked
                        )
                      }
                    />
                    Additional proteins
                  </label>

                  <Note>
                    STRING functional associations
                    are always included.
                    Overlapping links are drawn once.
                  </Note>
                </details>
              )}
            </>
          )}

          {/* WHOLE NETWORK */}

          {mode === "atlas" && (
            <>
              {networkFocus && (
                <div className="subset-panel">
                  <strong>
                    {networkFocus}: local network
                  </strong>

                  <Note>
                    Immediate neighbors within the selected
                    process and visible connection layers.
                  </Note>

                  <button
                    type="button"
                    className="connection-expand"
                    onClick={() => setNetworkFocus(null)}
                  >
                    Back to selected network
                  </button>
                </div>
              )}
              <Note>
                View the complete set of selected
                genes or isolate a biological process.
              </Note>

              <label className="field-label">
                Biological process
              </label>

              <select
                className="field"
                value={moduleId}
                onChange={e => {
                  setModuleId(e.target.value);
                  setNetworkFocus(null);
                }}
              >
                <option value="ALL">
                  All processes
                </option>

                {atlas.modules.map(m => (
                  <option
                    key={m.id}
                    value={m.id}
                  >
                    {m.label}
                  </option>
                ))}
              </select>

              <label className="mini-check">
                <input
                  type="checkbox"
                  checked={bridgesOnly}
                  onChange={e => {
                    setBridgesOnly(e.target.checked);
                    setNetworkFocus(null);
                  }}
                />

                Show connecting genes and their cell-fate links
              </label>

              <details className="layers">
                <summary>
                  Connection layers
                </summary>

                <label>
                  <input
                    type="checkbox"
                    checked={showDirect}
                    onChange={e =>
                      setShowDirect(
                        e.target.checked
                      )
                    }
                  />
                  Directed database interactions
                </label>

                <label>
                  <input
                    type="checkbox"
                    checked={showSaved}
                    onChange={e =>
                      setShowSaved(
                        e.target.checked
                      )
                    }
                  />
                  Earlier illustrative interactions
                </label>
              </details>

              <details className="compact-info">
                <summary>Interaction legend</summary>

                <div
                  style={{
                    display: "grid",
                    gap: "10px",
                    marginTop: "10px"
                  }}
                >
                  <div
                    style={{
                      display: "flex",
                      alignItems: "center",
                      gap: "10px"
                    }}
                  >
                    <span
                      style={{
                        display: "inline-block",
                        width: "28px",
                        height: "0",
                        borderTop: "3px solid rgba(180, 210, 255, 0.28)",
                        opacity: 0.95
                      }}
                    />
                    <span>STRING functional association</span>
                  </div>

                  <div
                    style={{
                      display: "flex",
                      alignItems: "center",
                      gap: "10px"
                    }}
                  >
                    <span
                      style={{
                        display: "inline-flex",
                        alignItems: "center",
                        width: "32px",
                        position: "relative"
                      }}
                    >
                      <span
                        style={{
                          display: "inline-block",
                          width: "24px",
                          height: "0",
                          borderTop: "3px solid #a8f1ff"
                        }}
                      />
                      <span
                        style={{
                          width: 0,
                          height: 0,
                          borderTop: "5px solid transparent",
                          borderBottom: "5px solid transparent",
                          borderLeft: "8px solid #a8f1ff"
                        }}
                      />
                    </span>
                    <span>Directed database interaction</span>
                  </div>

                  <div
                    style={{
                      display: "flex",
                      alignItems: "center",
                      gap: "10px"
                    }}
                  >
                    <span
                      style={{
                        display: "inline-flex",
                        alignItems: "center",
                        width: "32px",
                        position: "relative"
                      }}
                    >
                      <span
                        style={{
                          display: "inline-block",
                          width: "24px",
                          height: "0",
                          borderTop: "3px solid #ff9bac"
                        }}
                      />
                      <span
                        style={{
                          width: 0,
                          height: 0,
                          borderTop: "5px solid transparent",
                          borderBottom: "5px solid transparent",
                          borderLeft: "8px solid #ff9bac"
                        }}
                      />
                    </span>
                    <span>Inhibitory directed interaction</span>
                  </div>

                  <div
                    style={{
                      display: "flex",
                      alignItems: "center",
                      gap: "10px"
                    }}
                  >
                    <span
                      style={{
                        display: "inline-flex",
                        alignItems: "center",
                        width: "32px",
                        position: "relative"
                      }}
                    >
                      <span
                        style={{
                          display: "inline-block",
                          width: "24px",
                          height: "0",
                          borderTop: "3px solid #d4b3ff"
                        }}
                      />
                      <span
                        style={{
                          width: 0,
                          height: 0,
                          borderTop: "5px solid transparent",
                          borderBottom: "5px solid transparent",
                          borderLeft: "8px solid #d4b3ff"
                        }}
                      />
                    </span>
                    <span>Earlier illustrative interaction</span>
                  </div>
                </div>

                <p style={{ marginTop: "12px" }}>
                  Arrow direction comes from interaction databases.
                  It does not by itself prove that the same causal
                  regulation was active in the Fong neural stem-cell
                  experiment.
                </p>
              </details>
            </>
          )}

          {/* PATH EXPLORER */}

          {mode === "route" && (
            <>
              <div className="subset-switch">
                <button
                  className={
                    pathKind === "examples"
                      ? "active"
                      : ""
                  }
                  onClick={() => {
                    setPathKind("examples");
                    setRouteIndex(0);
                  }}
                >
                  Selected routes
                </button>

                <button
                  className={
                    pathKind === "saved"
                      ? "active"
                      : ""
                  }
                  onClick={() => {
                    setPathKind("saved");
                    setRouteIndex(0);
                  }}
                >
                  Search saved links
                </button>
              </div>

              {pathKind === "examples" ? (
                <>
                  <label className="field-label">
                    Start
                  </label>

                  <select
                    className="field"
                    value={source}
                    onChange={e => {
                      const next =
                        e.target.value;

                      setSource(next);

                      setTarget(
                        research.bypassRoutes.find(
                          r => r.source === next
                        )?.target || ""
                      );

                      setRouteIndex(0);
                    }}
                  >
                    {availableRouteSources.map(g => (
                      <option key={g}>
                        {g}
                      </option>
                    ))}
                  </select>

                  <label className="field-label">
                    Endpoint
                  </label>

                  <select
                    className="field"
                    value={target}
                    onChange={e => {
                      setTarget(e.target.value);
                      setRouteIndex(0);
                    }}
                  >
                    {availableRouteTargets.map(g => (
                      <option key={g}>
                        {g}
                      </option>
                    ))}
                  </select>

                  {activePath && (
                    <Note>
                      Example route avoiding{" "}
                      {activePath.blocker} in
                      the computer model.
                      This does not show a
                      mechanism operating
                      in neural cells.
                    </Note>
                  )}
                </>
              ) : (
                <>
                  <label className="field-label">
                    Start / endpoint
                  </label>


                  <div className="paired-fields">
                    <GeneAutocomplete
                      value={from}
                      names={geneNames}
                      placeholder="Starting gene"
                      onChange={value => {
                        setFrom(value);
                        setRouteIndex(0);
                      }}
                      onSelect={value => {
                        setFrom(value);
                        setRouteIndex(0);
                      }}
                    />

                    <GeneAutocomplete
                      value={to}
                      names={geneNames}
                      placeholder="Endpoint gene"
                      onChange={value => {
                        setTo(value);
                        setRouteIndex(0);
                      }}
                      onSelect={value => {
                        setTo(value);
                        setRouteIndex(0);
                      }}
                    />
                  </div>

                  <label className="field-label">
                    Path limits
                  </label>

                  <select
                    className="field"
                    value={setting}
                    onChange={e => {
                      setSetting(e.target.value);
                      setRouteIndex(0);
                    }}
                  >
                    {Object.entries(
                      atlas.metadata.settings
                    ).map(([id, rule]) => (
                      <option
                        key={id}
                        value={id}
                      >
                        Up to {rule.max_edges}
                        {" "}steps,{" "}
                        {rule.max_external}
                        {" "}extra proteins
                      </option>
                    ))}
                  </select>

                  <Note>
                    Searches only 393 saved example links,
                    not the complete interaction database.
                    A missing route is not evidence that
                    the genes are biologically disconnected.
                  </Note>
                </>
              )}

              {activePath ? (
                <>
                  <div className="route-stepper">
                    <button
                      disabled={!routeIndex}
                      onClick={() =>
                        setRouteIndex(i => i - 1)
                      }
                    >
                      ←
                    </button>

                    <span>
                      Route{" "}
                      {Math.min(
                        routeIndex + 1,
                        paths.length
                      )} of {paths.length}
                    </span>

                    <button
                      disabled={
                        routeIndex >=
                        paths.length - 1
                      }
                      onClick={() =>
                        setRouteIndex(i => i + 1)
                      }
                    >
                      →
                    </button>
                  </div>

                  <p className="route-note">
                    {pathNodes.join(" → ")}
                  </p>

                  {pathKind === "examples" && (
                    <details className="compact-info">
                      <summary>
                        Inspect each arrow
                      </summary>

                      {getRouteDetails.map(({
                        source: a,
                        target: b,
                        evidence,
                        warning
                      }) => (
                        <div
                          className="route-evidence"
                          key={edgeKey(a, b)}
                        >
                          <b>
                            {a} → {b}
                          </b>

                          <small>
                            {signText(
                              evidence?.effect
                            )}
                          </small>

                          {warning ? (
                            <p>{warning}</p>
                          ) : (
                            <p>
                              Biological direction
                              and cell-type relevance
                              not manually verified.
                            </p>
                          )}

                          {evidence?.references && (
                            <details>
                              <summary>
                                Reference identifiers
                              </summary>

                              <p>
                                {evidence.references}
                              </p>
                            </details>
                          )}
                        </div>
                      ))}
                    </details>
                  )}
                </>
              ) : (
                <Note>
                  No example available
                  for these settings.
                </Note>
              )}
            </>
          )}

        </aside>

        {/* RIGHT INSPECTOR */}

        {inspector && details && (
          <aside className="inspector compact-inspector">
            <button
              className="close"
              aria-label="Close details"
              onClick={() =>
                setInspector(false)
              }
            >
              ×
            </button>

            <div className="panel-kicker">
              SELECTED GENE
            </div>

            <h2>{details.id}</h2>

            <p
              className="kind"
              style={{
                color: getColor(details)
              }}
            >
              {details.tier?.includes("EXTERNAL")
                ? "Additional database protein"
                : details.f4_cross_module_bridge
                  ? "Connecting gene · related-process group"
                  : details.tier === "CORE"
                    ? "Cell-fate gene"
                    : "Related-process gene"}
            </p>

            {details.tier?.includes("EXTERNAL") && (
              <Note>
                This protein was added from prior-knowledge
                pathway records. It is not among the 321
                expression-selected genes; no RNA change
                is inferred from its presence here.
              </Note>
            )}

            {!details.tier?.includes("EXTERNAL") && (
              <>
                <div className="numbers">
                  <div>
                    <small>
                      RNA log₂ fold change
                    </small>

                    <strong>
                      {Number(
                        details.log2FC
                      ).toFixed(2)}
                    </strong>
                  </div>

                  <div>
                    <small>
                      Adjusted p-value
                    </small>

                    <strong>
                      {Number(
                        details.padj
                      ).toExponential(1)}
                    </strong>
                  </div>
                </div>

                <div className="tags">
                  {details.promoter_evidence && (
                    <span>
                      Promoter annotation
                    </span>
                  )}

                  {details.f4_cross_module_bridge && (
                    <span>
                      Connects groups
                    </span>
                  )}
                </div>

                <Note>
                  RNA is{" "}
                  {details.direction === "Up"
                    ? "higher"
                    : "lower"}{" "}
                  in the studied comparison.
                  This does not measure
                  protein activity.
                </Note>


                <h3 style={{ marginTop: 20, marginBottom: 10 }}>
                  Associated biological roles
                </h3>

                <div className="role-list">
                  {details.modules?.length ? (
                    details.modules.map(id => (
                      <div key={id}>
                        {modules.get(id) || id}
                      </div>
                    ))
                  ) : (
                    <div>No role annotation in this dataset</div>
                  )}
                </div>

                <p className="side-caption">
                  <b>{details.degree ?? "?"}</b> high-confidence
                  STRING functional associations in this atlas.
                  {details.tier === "INTERFACE" && (
                    <>
                      {" "}<b>{details.f3_core_neighbor_count ?? "?"}</b>
                      {" "}of these connect to the selected
                      cell-fate gene group.
                    </>
                  )}
                  {details.f4_cross_module_bridge && (
                    <> The gene connects at least two
                      network communities.</>
                  )}
                </p>

                {details.promoter_evidence && (
                  <p className="side-caption">
                    <b>Independent promoter evidence:</b>{" "}
                    this gene has a promoter annotation from
                    the Julien E2F binding comparison.
                    This does not prove direct regulation
                    in the Fong experiment.
                  </p>
                )}

                {research?.candidateIds?.includes(details.id) && (
                  <p className="side-caption">
                    <b>Closer investigation:</b>{" "}
                    included among the nine genes selected for
                    cross-study network and promoter review.
                  </p>
                )}

                {expression && (
                  <>
                    <h3 style={{ marginTop: 20, marginBottom: 10 }}>
                      Connectivity in the directed database
                    </h3>

                    {expression.mapped ? (
                      <div className="numbers">
                        <div>
                          <small>Reachable using selected genes</small>
                          <strong>{expression.outE4X0}</strong>
                        </div>

                        <div>
                          <small>Allowing one additional protein</small>
                          <strong>{expression.outE4X1}</strong>
                        </div>
                      </div>
                    ) : (
                      <Note>
                        No match in the directed database;
                        route analysis was not possible.
                      </Note>
                    )}
                  </>
                )}

                <details className="compact-info">
                  <summary>How to interpret this evidence</summary>
                  <p>
                    Biological roles and network communities are
                    annotations, not measurements of pathway
                    activation. STRING associations are
                    undirected. Reachability means that a
                    chain of at most four database arrows exists;
                    it does not establish regulation or cell death.
                  </p>
                </details>
              </>
            )}

            <button
              className="inspect-button"
              onClick={() => {
                pickFocus(details.id);
                setScope("choose");
              }}
            >
              Customize visible connections
            </button>

            <Note>
              Opens the neighbor picker in the left panel.
              Select only the genes you want to see with
              {" "}{details.id} in the 3D network.
            </Note>

            {mode === "focus" &&
              scope === "choose" &&
              focus === selected && (
                <div className="selected-brief">
                  <b>
                    {picked.length} selected partners
                  </b>

                  {picked.map(id => (
                    <button
                      key={id}
                      onClick={() =>
                        togglePick(id)
                      }
                    >
                      {id} ×
                    </button>
                  ))}
                </div>
              )}

            <details className="compact-info">
              <summary>
                Recorded interactions (
                {allSelectedEdges.length})
              </summary>

              <div className="inspector-neighbors">
                {[...new Set(
                  allSelectedEdges.map(e =>
                    e.source === selected
                      ? e.target
                      : e.source
                  )
                )]
                  .slice(0, 25)
                  .map(gene => (
                    <button
                      key={gene}
                      onClick={() =>
                        pickFocus(gene)
                      }
                    >
                      {gene}
                    </button>
                  ))}
              </div>

              <Note>
                These include undirected
                associations and published
                directional records;
                evidence strength varies.
              </Note>
            </details>
          </aside>
        )}

        {/* SCENE CONTROLS */}

        <div className="scene-tools">
          <button onClick={fit}>
            ◎ Fit
          </button>

          <button
            onClick={() =>
              graphRef.current?.cameraPosition(
                { x: 0, y: 0, z: 260 },
                { x: 0, y: 0, z: 0 },
                450
              )
            }
          >
            ↺ Reset view
          </button>
        </div>

        <div className="scene-hint">
          Rotate / Zoom / Click to inspect / Double-click for neighbors in this network
        </div>

        <div className="legend-bar">
          {[
            ["Cell-fate", COLORS.CORE],
            ["Related process", COLORS.INTERFACE],
            ["Connecting gene", COLORS.BRIDGE],
            ["Overlap annotation", COLORS.GOLD],
            ["Extra protein", COLORS.EXTERNAL]
          ].map(([label, color]) => (
            <span key={label}>
              <i
                style={{
                  background: color
                }}
              />
              {label}
            </span>
          ))}
        </div>
      </main>

      <footer>
        Expression data and network associations
        are distinct forms of evidence.
        Database routes are exploratory,
        not proven cell-death mechanisms.
      </footer>
    </div>
  );
}
