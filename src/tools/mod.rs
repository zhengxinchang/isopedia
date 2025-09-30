pub mod inspection;
pub mod manifest;
pub mod merge_replicates;

pub trait ToolCmdValidate {
    fn validate(&self) -> bool;
}
