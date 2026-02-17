pub enum Frame {
    First,
    Second,
    Third,
}

impl Frame {
    pub fn start_pos(&self) -> usize {
        match self {
            Frame::First => 0,
            Frame::Second => 1,
            Frame::Third => 2,
        }
    }
}
