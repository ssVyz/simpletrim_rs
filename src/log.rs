//! Console reporting on stderr, either as human readable text or as JSON lines.
//!
//! In JSON mode every line on stderr is one JSON object with an `event` field
//! (`start`, `info`, `adapter`, `progress`, `warning`, `error`, `summary`, ...) and the
//! elapsed time in seconds.

use std::io::Write;
use std::time::Instant;

use serde_json::{Map, Value, json};

pub struct Logger {
    json: bool,
    quiet: bool,
    start: Instant,
}

impl Logger {
    pub fn new(json: bool, quiet: bool) -> Self {
        Logger {
            json,
            quiet,
            start: Instant::now(),
        }
    }

    pub fn is_json(&self) -> bool {
        self.json
    }

    pub fn elapsed_secs(&self) -> f64 {
        self.start.elapsed().as_secs_f64()
    }

    fn write_line(&self, line: &str) {
        let mut err = std::io::stderr().lock();
        let _ = writeln!(err, "{line}");
    }

    fn emit_json(&self, event: &str, message: Option<&str>, fields: Value) {
        let mut obj = Map::new();
        obj.insert("event".into(), json!(event));
        obj.insert(
            "elapsed_s".into(),
            json!((self.elapsed_secs() * 1000.0).round() / 1000.0),
        );
        if let Some(m) = message {
            obj.insert("message".into(), json!(m));
        }
        if let Value::Object(extra) = fields {
            obj.extend(extra);
        }
        self.write_line(&Value::Object(obj).to_string());
    }

    /// Informational message; suppressed by `--quiet`.
    pub fn info(&self, message: &str) {
        self.info_with("info", message, Value::Null);
    }

    /// Informational message with a specific event name and extra JSON fields.
    pub fn info_with(&self, event: &str, message: &str, fields: Value) {
        if self.quiet {
            return;
        }
        if self.json {
            self.emit_json(event, Some(message), fields);
        } else {
            self.write_line(message);
        }
    }

    /// Progress update; suppressed by `--quiet`.
    pub fn progress(&self, message: &str, fields: Value) {
        self.info_with("progress", message, fields);
    }

    pub fn warn(&self, message: &str) {
        if self.json {
            self.emit_json("warning", Some(message), Value::Null);
        } else {
            self.write_line(&format!("warning: {message}"));
        }
    }

    pub fn error(&self, message: &str) {
        if self.json {
            self.emit_json("error", Some(message), Value::Null);
        } else {
            self.write_line(&format!("error: {message}"));
        }
    }

    /// Final result of a run; always printed. Human mode prints `text`, JSON mode emits
    /// `fields` as event `event`.
    pub fn result(&self, event: &str, text: &str, fields: Value) {
        if self.json {
            self.emit_json(event, None, fields);
        } else {
            self.write_line(text.trim_end());
        }
    }
}
