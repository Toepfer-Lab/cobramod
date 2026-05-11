// ── Interfaces ────────────────────────────────────────────────────────────────

export interface HoverData {
  display_name: string; id: string; kind_badge: string;
  comp_label: string; pipeline_tag: string;
  flux_str: string; std_str: string;
  substrates: string; products: string; extra: string;
}
export interface RxnNode {
  id: string; name: string; x: number; y: number;
  comp: string; kind: string; flux: number | null;
  fill_color: string; border_color: string;
  border_width: number; opacity: number; r: number;
  symbol: string; hover: HoverData;
}
export interface RailRouteView {
  line_id: string; width: number; color: string;
}
export interface StationLineView {
  klass: string;
  flux_sum: number;
  width: number;
  color: string;
}
export interface StationViewEntry {
  pair_id: string;
  lines: StationLineView[];
}
export interface ViewData {
  label: string; rxn_nodes: RxnNode[];
  rail_routes: RailRouteView[];
  stations: StationViewEntry[];
}
export interface RxnSummary {
  id: string; name: string; flux_per_view: string[];
  substrates?: string[]; products?: string[];
}
export interface RailRouteGeom {
  line_id: string; rxn_id: string; met_id: string;
  comp: string; stoich: number;
  base_color: string;
  seg_slots: number[];
  path: [number, number][];
}
export interface StationLineGeom {
  klass: string;
  rxn_ids: string[];
  rxn_summaries: RxnSummary[];
  spine_offset_idx: number;
  seg_slots: number[];
  path: [number, number][];
}
export interface StationGeom {
  pair_id: string;
  comp_a: string; comp_b: string;
  comp_a_label: string; comp_b_label: string;
  comp_a_color: string; comp_b_color: string;
  anchor_a: [number, number]; anchor_b: [number, number];
  tangent_a: [number, number]; tangent_b: [number, number];
  pill_height: number; pill_width: number;
  n_rxns: number;
  lines: StationLineGeom[];
}
export interface RouteDatum {
  pair_id: string;
  klass: string;
  rxn_ids: string[];
  seg_slots: number[];
  path: [number, number][];
  lineIndex: number;
}
export interface InnerEdge {
  rxn_id: string; met_id: string;
  x: [number, number]; y: [number, number];
}
export interface GuideLink {
  node_id: string; node_type: 'rxn' | 'met';
  x: [number, number]; y: [number, number];
}
export interface MetNode {
  id: string; name: string; x: number; y: number;
  comp: string; comp_label: string; comp_color: string;
  display_kind: 'stop' | 'detached';
  consumers: string[]; producers: string[];
}
export interface Compartment {
  key: string; label: string; color: string; fill: string;
  hull_vertices: [number, number][];
  hull_tension: number; label_x: number; label_y: number;
}
export interface Meta {
  model_name: string; x_range: [number, number]; y_range: [number, number];
  height: number; abs_max_flux: number; density_scale: number;
  grid_step: number;
}
export interface Interactivity {
  view_labels: string[];
  rxn_ids: string[]; rxn_names: string[]; rxn_x: number[]; rxn_y: number[];
  met_ids: string[]; met_names: string[]; met_x: number[]; met_y: number[];
}
export interface D3FigureData {
  meta: Meta; compartments: Compartment[];
  view_labels: string[]; views: ViewData[];
  met_nodes: MetNode[]; interactivity: Interactivity;
  rail_routes: RailRouteGeom[]; guide_links: GuideLink[];
  stations: StationGeom[];
  inner_edges: InnerEdge[];
}
export interface FocusState {
  type: 'rxn' | 'met' | 'station';
  id: string;
  rxnIds: Set<string>;
  metIds: Set<string>;
}

