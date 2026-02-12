//! Progress display utilities for CLI commands
//!
//! Provides honest progress indicators that show stages rather than
//! misleading percentage bars (which would freeze during large operations).

use indicatif::{ProgressBar, ProgressStyle};
use std::time::{Duration, Instant};

/// Memory usage tracker
#[allow(dead_code)]
pub struct MemoryTracker {
    start_memory: usize,
    peak_memory: usize,
}

impl MemoryTracker {
    pub fn new() -> Self {
        let current = Self::get_current_memory();
        Self {
            start_memory: current,
            peak_memory: current,
        }
    }

    /// Update peak memory if current is higher
    pub fn update(&mut self) {
        let current = Self::get_current_memory();
        if current > self.peak_memory {
            self.peak_memory = current;
        }
    }

    /// Get current memory usage in bytes (approximate)
    fn get_current_memory() -> usize {
        // Use a simple heuristic based on allocated memory
        // This is platform-dependent and approximate
        #[cfg(target_os = "linux")]
        {
            if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
                for line in status.lines() {
                    if line.starts_with("VmRSS:") {
                        if let Some(kb_str) = line.split_whitespace().nth(1) {
                            if let Ok(kb) = kb_str.parse::<usize>() {
                                return kb * 1024;
                            }
                        }
                    }
                }
            }
        }
        // Fallback: return 0 if we can't determine memory
        0
    }

    /// Get peak memory usage in bytes
    #[allow(dead_code)]
    pub fn peak_bytes(&self) -> usize {
        self.peak_memory
    }

    /// Get peak memory usage formatted as string
    pub fn peak_formatted(&self) -> String {
        format_bytes(self.peak_memory)
    }

    /// Get memory increase since start
    #[allow(dead_code)]
    pub fn increase_formatted(&self) -> String {
        let increase = self.peak_memory.saturating_sub(self.start_memory);
        format_bytes(increase)
    }
}

impl Default for MemoryTracker {
    fn default() -> Self {
        Self::new()
    }
}

/// Format bytes as human-readable string
pub fn format_bytes(bytes: usize) -> String {
    const KB: usize = 1024;
    const MB: usize = KB * 1024;
    const GB: usize = MB * 1024;

    if bytes >= GB {
        format!("{:.2} GB", bytes as f64 / GB as f64)
    } else if bytes >= MB {
        format!("{:.2} MB", bytes as f64 / MB as f64)
    } else if bytes >= KB {
        format!("{:.2} KB", bytes as f64 / KB as f64)
    } else {
        format!("{} B", bytes)
    }
}

/// Format duration as human-readable string
pub fn format_duration(duration: Duration) -> String {
    let secs = duration.as_secs();
    if secs >= 3600 {
        let hours = secs / 3600;
        let mins = (secs % 3600) / 60;
        format!("{}h {}m", hours, mins)
    } else if secs >= 60 {
        let mins = secs / 60;
        let secs = secs % 60;
        format!("{}m {}s", mins, secs)
    } else if secs > 0 {
        format!("{}s", secs)
    } else {
        format!("{}ms", duration.as_millis())
    }
}

/// Stage-based progress tracker
/// 
/// Shows honest progress by displaying stages rather than misleading percentages.
/// This is more appropriate for operations where we can't predict exact progress.
pub struct StageProgress {
    total_stages: usize,
    current_stage: usize,
    stage_names: Vec<String>,
    stage_times: Vec<Duration>,
    start_time: Instant,
    stage_start: Instant,
    spinner: ProgressBar,
    memory_tracker: MemoryTracker,
    verbose: bool,
}

impl StageProgress {
    /// Create a new stage progress tracker
    pub fn new(stages: Vec<&str>, verbose: bool) -> Self {
        let total_stages = stages.len();
        let stage_names: Vec<String> = stages.into_iter().map(|s| s.to_string()).collect();
        
        let spinner = ProgressBar::new_spinner();
        spinner.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.cyan} {msg}")
                .unwrap(),
        );
        spinner.enable_steady_tick(Duration::from_millis(100));

        Self {
            total_stages,
            current_stage: 0,
            stage_names,
            stage_times: Vec::new(),
            start_time: Instant::now(),
            stage_start: Instant::now(),
            spinner,
            memory_tracker: MemoryTracker::new(),
            verbose,
        }
    }

    /// Start the next stage
    pub fn next_stage(&mut self) {
        // Record time for previous stage
        if self.current_stage > 0 {
            let elapsed = self.stage_start.elapsed();
            self.stage_times.push(elapsed);
        }

        self.current_stage += 1;
        self.stage_start = Instant::now();
        self.memory_tracker.update();

        if self.current_stage <= self.total_stages {
            let stage_name = &self.stage_names[self.current_stage - 1];
            let msg = format!(
                "[{}/{}] ⏳ {}",
                self.current_stage, self.total_stages, stage_name
            );
            self.spinner.set_message(msg);
        }
    }

    /// Mark current stage as complete with optional details
    pub fn complete_stage(&mut self, details: Option<&str>) {
        self.memory_tracker.update();
        let elapsed = self.stage_start.elapsed();
        
        if self.current_stage <= self.total_stages {
            let stage_name = &self.stage_names[self.current_stage - 1];
            let time_str = format_duration(elapsed);
            
            let msg = if let Some(detail) = details {
                format!(
                    "[{}/{}] ✓ {} ({}) [{}]",
                    self.current_stage, self.total_stages, stage_name, detail, time_str
                )
            } else {
                format!(
                    "[{}/{}] ✓ {} [{}]",
                    self.current_stage, self.total_stages, stage_name, time_str
                )
            };
            
            self.spinner.suspend(|| {
                println!("{}", msg);
            });
        }
    }

    /// Update the current stage message (for sub-progress)
    pub fn update_message(&self, sub_message: &str) {
        if self.current_stage <= self.total_stages {
            let stage_name = &self.stage_names[self.current_stage - 1];
            let elapsed = self.stage_start.elapsed();
            let time_str = format_duration(elapsed);
            
            let msg = format!(
                "[{}/{}] ⏳ {} - {} [{}]",
                self.current_stage, self.total_stages, stage_name, sub_message, time_str
            );
            self.spinner.set_message(msg);
        }
    }

    /// Print a message without interrupting the spinner
    pub fn println(&self, message: &str) {
        self.spinner.suspend(|| {
            println!("{}", message);
        });
    }

    /// Print verbose message (only if verbose mode is enabled)
    pub fn verbose(&self, message: &str) {
        if self.verbose {
            self.spinner.suspend(|| {
                println!("   {}", message);
            });
        }
    }

    /// Finish all stages and show summary
    pub fn finish(&mut self) {
        // Record final stage time
        if self.current_stage > 0 && self.stage_times.len() < self.current_stage {
            let elapsed = self.stage_start.elapsed();
            self.stage_times.push(elapsed);
        }

        self.memory_tracker.update();
        self.spinner.finish_and_clear();

        let total_time = self.start_time.elapsed();
        let peak_mem = self.memory_tracker.peak_formatted();

        println!();
        println!("✅ Completed in {} | Peak memory: {}", 
            format_duration(total_time), peak_mem);
    }

    /// Finish with error
    #[allow(dead_code)]
    pub fn finish_with_error(&mut self, error: &str) {
        self.spinner.finish_and_clear();
        
        let total_time = self.start_time.elapsed();
        println!();
        println!("❌ Failed at stage {}/{}: {}", 
            self.current_stage, self.total_stages, error);
        println!("   Elapsed time: {}", format_duration(total_time));
    }

    /// Get total elapsed time
    #[allow(dead_code)]
    pub fn elapsed(&self) -> Duration {
        self.start_time.elapsed()
    }

    /// Get peak memory usage
    #[allow(dead_code)]
    pub fn peak_memory(&self) -> String {
        self.memory_tracker.peak_formatted()
    }
}

/// Simple spinner for operations without stages
#[allow(dead_code)]
pub struct SimpleSpinner {
    spinner: ProgressBar,
    start_time: Instant,
    memory_tracker: MemoryTracker,
}

#[allow(dead_code)]
impl SimpleSpinner {
    pub fn new(message: &str) -> Self {
        let spinner = ProgressBar::new_spinner();
        spinner.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.cyan} {msg}")
                .unwrap(),
        );
        spinner.set_message(message.to_string());
        spinner.enable_steady_tick(Duration::from_millis(100));

        Self {
            spinner,
            start_time: Instant::now(),
            memory_tracker: MemoryTracker::new(),
        }
    }

    pub fn set_message(&self, message: &str) {
        let elapsed = format_duration(self.start_time.elapsed());
        self.spinner.set_message(format!("{} [{}]", message, elapsed));
    }

    pub fn finish(&mut self, message: &str) {
        self.memory_tracker.update();
        self.spinner.finish_and_clear();
        
        let elapsed = format_duration(self.start_time.elapsed());
        let peak_mem = self.memory_tracker.peak_formatted();
        
        println!("{} [{}] | Peak memory: {}", message, elapsed, peak_mem);
    }

    pub fn finish_with_error(&mut self, error: &str) {
        self.spinner.finish_and_clear();
        println!("❌ {}", error);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_bytes() {
        assert_eq!(format_bytes(500), "500 B");
        assert_eq!(format_bytes(1024), "1.00 KB");
        assert_eq!(format_bytes(1024 * 1024), "1.00 MB");
        assert_eq!(format_bytes(1024 * 1024 * 1024), "1.00 GB");
        assert_eq!(format_bytes(1536 * 1024 * 1024), "1.50 GB");
    }

    #[test]
    fn test_format_duration() {
        assert_eq!(format_duration(Duration::from_millis(500)), "500ms");
        assert_eq!(format_duration(Duration::from_secs(30)), "30s");
        assert_eq!(format_duration(Duration::from_secs(90)), "1m 30s");
        assert_eq!(format_duration(Duration::from_secs(3700)), "1h 1m");
    }
}
