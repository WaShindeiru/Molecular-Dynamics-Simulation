use std::path::{Path, PathBuf};

/// Returns `dir/file_name` unless that path already exists, in which case
/// returns `dir/{stem}_end.json` (e.g. `parameters.json` → `parameters_end.json`).
pub fn json_save_path_avoiding_overwrite(dir: &Path, file_name: &str) -> PathBuf {
  let path = dir.join(file_name);

  if path.exists() {
    let stem = file_name.strip_suffix(".json").unwrap_or(file_name);
    dir.join(format!("{stem}_end.json"))
  } else {
    path
  }
}

#[cfg(test)]
mod tests {
  use std::fs;

  use super::*;

  #[test]
  fn uses_default_name_when_file_missing() {
    let dir = std::env::temp_dir().join(format!(
      "json_save_path_test_{}",
      std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();

    let path = json_save_path_avoiding_overwrite(&dir, "parameters.json");
    assert_eq!(path, dir.join("parameters.json"));

    let _ = fs::remove_dir_all(&dir);
  }

  #[test]
  fn uses_end_suffix_when_file_exists() {
    let dir = std::env::temp_dir().join(format!(
      "json_save_path_exists_test_{}",
      std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();

    let existing = dir.join("generator_config.json");
    fs::write(&existing, "{}").unwrap();

    let path = json_save_path_avoiding_overwrite(&dir, "generator_config.json");
    assert_eq!(path, dir.join("generator_config_end.json"));

    let _ = fs::remove_dir_all(&dir);
  }
}
