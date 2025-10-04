pub mod manifest;
pub mod output;
pub mod profile;
pub mod view;

pub trait ToolCmdValidate {
    fn validate(&self) -> bool;
}
