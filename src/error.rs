use std::fmt;
use std::io;
use std::num::ParseIntError;
#[derive(Debug)]
pub enum MyError {
    Io(io::Error),
    ParseInt(ParseIntError),
    InvalidInput(String),
    Serde(String),
    EmptyRes(String),
}

impl fmt::Display for MyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MyError::Io(err) => write!(f, "IO error: {}", err),
            MyError::ParseInt(err) => write!(f, "Parse int error: {}", err),
            MyError::InvalidInput(msg) => write!(f, "InvalidInput: {}", msg),
            MyError::Serde(msg) => write!(f, "Serde error: {}", msg),
            MyError::EmptyRes(msg) => write!(f, "Empty result: {}", msg),
        }
    }
}

impl std::error::Error for MyError {}

impl From<io::Error> for MyError {
    fn from(err: io::Error) -> MyError {
        MyError::Io(err)
    }
}

impl From<ParseIntError> for MyError {
    fn from(err: ParseIntError) -> MyError {
        MyError::ParseInt(err)
    }
}

impl From<serde_json::Error> for MyError {
    fn from(err: serde_json::Error) -> MyError {
        MyError::Serde(err.to_string())
    }
}
