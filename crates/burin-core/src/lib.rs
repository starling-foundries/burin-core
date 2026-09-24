//! A deterministic, equal-area coverage commitment: polygon → rHEALPix cells → canonical cell
//! set → 32-byte root, with openings (membership / non-membership), set algebra over roots, and a
//! profile that binds every parameter the identifier is read against; and the zone topology
//! (scanline sub-zone order, edge neighbours, rasters) that never changes a root.

pub mod error;
pub mod geo;
pub mod hash;
pub mod hierarchy;
pub mod index;
pub mod opening;
pub mod polyfill;
pub mod profile;
pub mod setops;
pub mod tree;
pub mod zone;
pub mod zone_data;

pub use error::{Error, Result};
pub use hash::{Ctx, Digest, Hasher, Ladders, Sha256Tagged};
pub use hierarchy::{cid_to_suid, suid_to_cid, Cid, Hierarchy, MAX_APERTURE, SPACE};
pub use index::{Hit, Index, Relation};
pub use opening::{open_path, verify_opening, Claim, Entry, Opening, OpeningRecord};
pub use polyfill::{cover, cover_nodes, fingerprint_geojson, fingerprint_polygon, geometry_from_geojson, polyfill, MAX_RESOLUTION};
pub use profile::{Profile, BURIN_1, OGC_RHEALPIX};
pub use setops::{difference, divergence, intersect, merge, prove, union, Op, SetOpProof};
pub use tree::{Node, Tree};
pub use zone_data::{dggrs_uri, from_dggs_json, to_dggs_json, Presence, RHEALPIX_DGGRS};
pub use zone::{cell_at, halo_index, neighbours, position, raster_cells, raster_index, subzone_at, subzone_index, subzones, Direction, Position, DIRECTIONS};
