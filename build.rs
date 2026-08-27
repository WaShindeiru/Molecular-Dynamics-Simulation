use std::process::Command;

fn main() {
  // Rebuild when HEAD or the index changes (covers new commits and dirty tree).
  println!("cargo:rerun-if-changed=.git/HEAD");
  println!("cargo:rerun-if-changed=.git/index");

  let hash = git_output(&["rev-parse", "HEAD"]).unwrap_or_else(|| "unknown".to_string());
  let dirty = git_output(&["status", "--porcelain"])
    .map(|s| !s.is_empty())
    .unwrap_or(false);

  let git_hash = if dirty && hash != "unknown" {
    format!("{hash}-dirty")
  } else {
    hash
  };

  println!("cargo:rustc-env=GIT_HASH={git_hash}");
}

fn git_output(args: &[&str]) -> Option<String> {
  let output = Command::new("git").args(args).output().ok()?;
  if !output.status.success() {
    return None;
  }
  let s = String::from_utf8(output.stdout).ok()?;
  Some(s.trim().to_string())
}
