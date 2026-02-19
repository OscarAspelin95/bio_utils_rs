use std::collections::HashSet;

pub trait Sketcher: Send + Sync {
    fn sketch(&self, seq: &[u8]) -> HashSet<u64>;
}
