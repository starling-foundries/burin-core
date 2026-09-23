//! The hierarchy H(A, B) and the cell id.
//!
//! ```text
//! path   (b, d_1, ..., d_r)     b in [0, B), d_i in [0, A), r >= 0 the level
//! cid    (A + b) * A^r + Σ d_i * A^(r-1-i)     one integer
//! suid   'Q453'                                the rHEALPix human form, B == 6 only
//! ```
//! With `B <= A(A-1)` a cid at level r has exactly r + 2 base-A digits, so `level`, `parent`
//! (`cid / A`), `child` (`A*cid + k`) and contiguous `descendants` all follow from the digits.

use crate::error::{invalid, Result};

pub type Cid = u64;

pub const BASE_LETTERS: &[u8; 6] = b"NOPQRS";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Hierarchy {
    pub a: u32,
    pub b: u32,
}

/// rHEALPix with `N_side = 3`: aperture 9 over six base cells.
pub const SPACE: Hierarchy = Hierarchy { a: 9, b: 6 };

impl Hierarchy {
    pub fn new(a: u32, b: u32) -> Result<Hierarchy> {
        if a < 2 || b < 1 {
            return invalid(format!("need A >= 2 and B >= 1, got A={a}, B={b}"));
        }
        if b > a * (a - 1) {
            return invalid(format!("B={b} > A(A-1)={}: the leading block would not be two digits", a * (a - 1)));
        }
        Ok(Hierarchy { a, b })
    }

    pub fn cid(&self, path: &[u32]) -> Result<Cid> {
        let (&b, digits) = match path.split_first() {
            Some(p) => p,
            None => return invalid("empty path"),
        };
        if b >= self.b {
            return invalid(format!("base {b} out of range [0, {})", self.b));
        }
        let mut c = (self.a + b) as u64;
        for &d in digits {
            if d >= self.a {
                return invalid(format!("digit {d} out of range [0, {})", self.a));
            }
            c = c * self.a as u64 + d as u64;
        }
        Ok(c)
    }

    pub fn level(&self, cid: Cid) -> u32 {
        let (mut n, mut x) = (0u32, cid);
        while x != 0 {
            x /= self.a as u64;
            n += 1;
        }
        n.saturating_sub(2)
    }

    pub fn path(&self, cid: Cid) -> Result<Vec<u32>> {
        let r = self.check(cid, None)?;
        let mut digits = Vec::with_capacity(r as usize + 1);
        let mut x = cid;
        for _ in 0..r {
            digits.push((x % self.a as u64) as u32);
            x /= self.a as u64;
        }
        digits.push((x - self.a as u64) as u32);
        digits.reverse();
        Ok(digits)
    }

    pub fn parent(&self, cid: Cid) -> Option<Cid> {
        if self.level(cid) == 0 {
            None
        } else {
            Some(cid / self.a as u64)
        }
    }

    pub fn child(&self, cid: Cid, k: u32) -> Cid {
        cid * self.a as u64 + k as u64
    }

    pub fn base(&self, cid: Cid) -> u32 {
        (cid / (self.a as u64).pow(self.level(cid)) - self.a as u64) as u32
    }

    /// True iff `a` is `b` or a proper ancestor of `b`.
    pub fn is_ancestor(&self, a: Cid, b: Cid) -> bool {
        let (la, lb) = (self.level(a), self.level(b));
        la <= lb && b / (self.a as u64).pow(lb - la) == a
    }

    /// The contiguous cid range of the descendants `k` levels below `cid`.
    pub fn descendants(&self, cid: Cid, k: u32) -> core::ops::Range<Cid> {
        let span = (self.a as u64).pow(k);
        cid * span..(cid + 1) * span
    }

    /// Validate a cid (well-formed leading block, level bound); returns its level.
    pub fn check(&self, cid: Cid, max_level: Option<u32>) -> Result<u32> {
        if cid < self.a as u64 {
            return invalid(format!("not a cid: {cid}"));
        }
        let r = self.level(cid);
        let lead = cid / (self.a as u64).pow(r);
        if !(self.a as u64 <= lead && lead < (self.a + self.b) as u64) {
            return invalid(format!("cid {cid} names base {}, outside [0, {})", lead as i64 - self.a as i64, self.b));
        }
        if let Some(m) = max_level {
            if r > m {
                return invalid(format!("cid {cid} is at level {r}, deeper than {m}"));
            }
        }
        Ok(r)
    }
}

/// `'Q453'` → `[3, 4, 5, 3]`.
pub fn suid_to_path(suid: &str) -> Result<Vec<u32>> {
    let bytes = suid.as_bytes();
    let (&head, rest) = match bytes.split_first() {
        Some(p) => p,
        None => return invalid("empty suid"),
    };
    let b = match BASE_LETTERS.iter().position(|&c| c == head) {
        Some(b) => b as u32,
        None => return invalid(format!("unknown base cell {:?}", head as char)),
    };
    let mut path = vec![b];
    for &c in rest {
        if !c.is_ascii_digit() {
            return invalid(format!("bad suid digit {:?}", c as char));
        }
        path.push((c - b'0') as u32);
    }
    Ok(path)
}

/// `[3, 4, 5, 3]` → `'Q453'`.
pub fn path_to_suid(path: &[u32]) -> String {
    let mut s = String::with_capacity(path.len());
    s.push(BASE_LETTERS[path[0] as usize] as char);
    for &d in &path[1..] {
        s.push((b'0' + d as u8) as char);
    }
    s
}

pub fn suid_to_cid(suid: &str) -> Result<Cid> {
    SPACE.cid(&suid_to_path(suid)?)
}

pub fn cid_to_suid(cid: Cid) -> Result<String> {
    Ok(path_to_suid(&SPACE.path(cid)?))
}
