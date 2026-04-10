import csv
import json
import queue
import threading
import traceback
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

from pipeline_ABCDES import pipeline_ABCDES


APP_TITLE = "PrimerScore"
DEFAULT_ORGANISM = "Glycine max"
DEFAULT_TOP_N = 10
DEFAULT_DUMMY = True

DEMO_SEQUENCE = """>demo_sequence
ATGGAGAAAGTTGGAATGTTGTTGCTGCTGCTGCTGCTGCTGGTGCTGATGATGATGCTGCT
GCTGATGCTGCTGATGATGCTGCTGATGCTGATGCTGCTGATGCTGCTGATGATGCTGATGC
TGCTGATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGATGCTGATGCTGCTGATGCTG
ATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGATGCT
GATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGATGC
TGATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGATG
CTGATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGAT
GCTGATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTGATGCTGATGCTGATGCTGCTTAA
"""


class PrimerScoreGUI:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title(APP_TITLE)
        self.root.geometry("1320x860")
        self.root.minsize(1120, 760)

        self.status_var = tk.StringVar(value="Ready")
        self.organism_var = tk.StringVar(value=DEFAULT_ORGANISM)
        self.top_n_var = tk.StringVar(value=str(DEFAULT_TOP_N))
        self.dummy_var = tk.BooleanVar(value=DEFAULT_DUMMY)

        self.result_data = None
        self.current_rows = []
        self.worker_thread = None
        self.log_queue = queue.Queue()

        self._build_layout()
        self._poll_log_queue()

    # -------------------------
    # UI Layout
    # -------------------------
    def _build_layout(self):
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(3, weight=1)

        self._build_header()
        self._build_top_panels()
        self._build_log_panel()
        self._build_result_panel()

    def _build_header(self):
        frame = ttk.Frame(self.root, padding=(12, 10))
        frame.grid(row=0, column=0, sticky="ew")
        frame.columnconfigure(0, weight=1)

        ttk.Label(frame, text="PrimerScore", font=("Arial", 18, "bold")).grid(
            row=0, column=0, sticky="w"
        )
        ttk.Label(
            frame,
            text="PCR primer design with BLAST, expression-aware evaluation and scoring",
            font=("Arial", 10),
        ).grid(row=1, column=0, sticky="w", pady=(2, 0))

        ttk.Label(
            frame,
            textvariable=self.status_var,
            font=("Arial", 10, "bold"),
            foreground="#1f4e79",
        ).grid(row=0, column=1, rowspan=2, sticky="e")

    def _build_top_panels(self):
        top = ttk.Frame(self.root, padding=(12, 0, 12, 8))
        top.grid(row=1, column=0, sticky="nsew")
        top.columnconfigure(0, weight=3)
        top.columnconfigure(1, weight=1)

        self._build_input_panel(top)
        self._build_control_panel(top)

    def _build_input_panel(self, parent):
        frame = ttk.LabelFrame(parent, text="Target Sequence", padding=10)
        frame.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(1, weight=1)

        ttk.Label(
            frame,
            text="Paste raw DNA sequence or FASTA. Non-ACGT characters will be removed automatically.",
        ).grid(row=0, column=0, sticky="w", pady=(0, 6))

        self.sequence_text = ScrolledText(
            frame, wrap=tk.WORD, height=12, font=("Consolas", 10)
        )
        self.sequence_text.grid(row=1, column=0, sticky="nsew")

    def _build_control_panel(self, parent):
        frame = ttk.LabelFrame(parent, text="Run Settings", padding=10)
        frame.grid(row=0, column=1, sticky="nsew")
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="Organism").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Entry(frame, textvariable=self.organism_var).grid(
            row=0, column=1, sticky="ew", pady=4
        )

        ttk.Label(frame, text="Top N").grid(row=1, column=0, sticky="w", pady=4)
        ttk.Entry(frame, textvariable=self.top_n_var).grid(
            row=1, column=1, sticky="ew", pady=4
        )

        ttk.Checkbutton(
            frame, text="Use dummy BLAST mode", variable=self.dummy_var
        ).grid(row=2, column=0, columnspan=2, sticky="w", pady=6)

        ttk.Separator(frame, orient="horizontal").grid(
            row=3, column=0, columnspan=2, sticky="ew", pady=8
        )

        self.run_button = ttk.Button(frame, text="Run PrimerScore", command=self._on_run)
        self.run_button.grid(row=4, column=0, columnspan=2, sticky="ew", pady=(2, 6))

        ttk.Button(frame, text="Load Demo", command=self._on_load_demo).grid(
            row=5, column=0, columnspan=2, sticky="ew", pady=3
        )
        ttk.Button(frame, text="Clear", command=self._on_clear).grid(
            row=6, column=0, columnspan=2, sticky="ew", pady=3
        )

        ttk.Label(
            frame,
            text="Tip: first test with dummy mode.\nAfter GUI works, uncheck it for real BLAST.",
            foreground="#555555",
            justify="left",
        ).grid(row=7, column=0, columnspan=2, sticky="w", pady=(10, 0))

    def _build_log_panel(self):
        frame = ttk.LabelFrame(self.root, text="Log Output", padding=10)
        frame.grid(row=2, column=0, sticky="nsew", padx=12, pady=(0, 8))
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        self.log_text = ScrolledText(
            frame, wrap=tk.WORD, height=10, state="disabled", font=("Consolas", 10)
        )
        self.log_text.grid(row=0, column=0, sticky="nsew")

    def _build_result_panel(self):
        frame = ttk.LabelFrame(self.root, text="Combined Ranking Results", padding=10)
        frame.grid(row=3, column=0, sticky="nsew", padx=12, pady=(0, 12))
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        columns = (
            "rank",
            "pair_index",
            "forward",
            "reverse",
            "product_size",
            "source",
            "score",
            "label",
            "warnings",
        )

        self.tree = ttk.Treeview(frame, columns=columns, show="headings", height=18)
        self.tree.grid(row=0, column=0, sticky="nsew")
        self.tree.bind("<Double-1>", self._on_row_double_click)

        widths = {
            "rank": 60,
            "pair_index": 90,
            "forward": 240,
            "reverse": 240,
            "product_size": 100,
            "source": 120,
            "score": 90,
            "label": 120,
            "warnings": 260,
        }

        for col in columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=widths.get(col, 100), anchor="w")

        y_scroll = ttk.Scrollbar(frame, orient="vertical", command=self.tree.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        self.tree.configure(yscrollcommand=y_scroll.set)

        bottom = ttk.Frame(frame)
        bottom.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        bottom.columnconfigure(0, weight=1)

        self.result_hint = ttk.Label(bottom, text="No results yet.")
        self.result_hint.grid(row=0, column=0, sticky="w")

        right = ttk.Frame(bottom)
        right.grid(row=0, column=1, sticky="e")

        self.export_json_button = ttk.Button(
            right, text="Export JSON", command=self._on_export_json, state="disabled"
        )
        self.export_json_button.grid(row=0, column=0, padx=4)

        self.export_csv_button = ttk.Button(
            right, text="Export CSV", command=self._on_export_csv, state="disabled"
        )
        self.export_csv_button.grid(row=0, column=1, padx=4)

    # -------------------------
    # Logging / Status
    # -------------------------
    def _append_log(self, message: str):
        self.log_text.configure(state="normal")
        self.log_text.insert(tk.END, message.rstrip() + "\n")
        self.log_text.see(tk.END)
        self.log_text.configure(state="disabled")

    def _queue_log(self, message: str):
        self.log_queue.put(("log", message))

    def _queue_status(self, status: str):
        self.log_queue.put(("status", status))

    def _poll_log_queue(self):
        try:
            while True:
                item_type, value = self.log_queue.get_nowait()
                if item_type == "log":
                    self._append_log(value)
                elif item_type == "status":
                    self.status_var.set(value)
        except queue.Empty:
            pass
        self.root.after(100, self._poll_log_queue)

    # -------------------------
    # Input / Clear / Demo
    # -------------------------
    def _on_load_demo(self):
        self.sequence_text.delete("1.0", tk.END)
        self.sequence_text.insert("1.0", DEMO_SEQUENCE)
        self._append_log("[INFO] Demo sequence loaded.")

    def _on_clear(self):
        self.sequence_text.delete("1.0", tk.END)
        self._clear_results()
        self.result_data = None
        self.current_rows = []
        self.status_var.set("Ready")

        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", tk.END)
        self.log_text.configure(state="disabled")

    def _clear_results(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.result_hint.config(text="No results yet.")
        self.export_json_button.config(state="disabled")
        self.export_csv_button.config(state="disabled")

    def _collect_inputs(self):
        raw_sequence = self.sequence_text.get("1.0", tk.END).strip()
        sequence = self._clean_sequence(raw_sequence)
        if not sequence:
            raise ValueError("Sequence is empty after cleaning.")
        if len(sequence) < 50:
            raise ValueError("Sequence is too short. Please provide a longer target sequence.")

        organism = self.organism_var.get().strip() or DEFAULT_ORGANISM

        top_n_text = self.top_n_var.get().strip()
        if not top_n_text:
            raise ValueError("Top N cannot be empty.")
        try:
            top_n = int(top_n_text)
        except ValueError as exc:
            raise ValueError("Top N must be an integer.") from exc
        if top_n <= 0:
            raise ValueError("Top N must be greater than 0.")

        return {
            "sequence": sequence,
            "organism": organism,
            "top_n": top_n,
            "dummy_blast": self.dummy_var.get(),
        }

    @staticmethod
    def _clean_sequence(raw_text: str) -> str:
        lines = [line.strip() for line in raw_text.splitlines() if line.strip()]
        if not lines:
            return ""
        if lines[0].startswith(">"):
            lines = lines[1:]
        joined = "".join(lines).upper()
        return "".join(ch for ch in joined if ch in {"A", "C", "G", "T"})

    # -------------------------
    # Run
    # -------------------------
    def _on_run(self):
        if self.worker_thread and self.worker_thread.is_alive():
            messagebox.showinfo(APP_TITLE, "A run is already in progress.")
            return

        try:
            inputs = self._collect_inputs()
        except Exception as exc:
            messagebox.showerror(APP_TITLE, str(exc))
            return

        self._clear_results()
        self.result_data = None
        self.current_rows = []

        self.run_button.config(state="disabled")
        self._queue_status("Running")
        self._queue_log("[INFO] PrimerScore started.")
        self._queue_log(f"[INFO] Sequence length: {len(inputs['sequence'])}")
        self._queue_log(f"[INFO] Organism: {inputs['organism']}")
        self._queue_log(f"[INFO] Dummy BLAST: {inputs['dummy_blast']}")
        self._queue_log("[INFO] Launching pipeline thread...")

        self.worker_thread = threading.Thread(
            target=self._run_pipeline_worker, args=(inputs,), daemon=True
        )
        self.worker_thread.start()

    def _run_pipeline_worker(self, inputs: dict):
        try:
            self._queue_log("[INFO] Calling pipeline_ABCDES(...)")
            result = pipeline_ABCDES(
                sequence=inputs["sequence"],
                organism=inputs["organism"],
                dummy_blast=inputs["dummy_blast"],
                enable_e=True,
                enable_scoring=True,
            )

            meta = result.get("meta", {})
            counts = meta.get("counts", {})
            e_status = meta.get("e_status", {})
            scoring_status = meta.get("scoring_status", {})

            self._queue_log(
                f"[INFO] Counts: pairs={counts.get('pairs', 0)}, "
                f"blast_results={counts.get('blast_results', 0)}, "
                f"combined_ranking={counts.get('combined_ranking', 0)}"
            )
            self._queue_log(
                f"[INFO] E status: ok={e_status.get('ok')}, "
                f"resolved_target_pairs={e_status.get('resolved_target_pairs', 0)}, "
                f"pairs_with_nonempty_target_expr={e_status.get('pairs_with_nonempty_target_expr', 0)}"
            )
            self._queue_log(
                f"[INFO] Scoring status: ok={scoring_status.get('ok')}, "
                f"n_pairs_in_d={scoring_status.get('n_pairs_in_d', 0)}, "
                f"n_pairs_in_e={scoring_status.get('n_pairs_in_e', 0)}"
            )

            rows = self._ranking_to_rows(result, inputs["top_n"])
            self.root.after(0, lambda: self._update_ui_success(result, rows))
        except Exception:
            err = traceback.format_exc()
            self.root.after(0, lambda: self._update_ui_failure(err))

    def _update_ui_success(self, result: dict, rows: list):
        self.result_data = result
        self.current_rows = rows
        self._fill_result_table(rows)

        self.run_button.config(state="normal")
        self.status_var.set("Finished")
        self._append_log("[INFO] Pipeline finished successfully.")

        if self.result_data:
            self.export_json_button.config(state="normal")
        if self.current_rows:
            self.export_csv_button.config(state="normal")

    def _update_ui_failure(self, error_text: str):
        self.run_button.config(state="normal")
        self.status_var.set("Failed")
        self._append_log("[ERROR] Pipeline failed.")
        self._append_log(error_text)
        messagebox.showerror(APP_TITLE, "Run failed. See log for details.")

    # -------------------------
    # Data lookup helpers
    # -------------------------
    @staticmethod
    def _normalize_pair_index(value):
        try:
            return int(value)
        except Exception:
            return value

    @staticmethod
    def _extract_sequence_from_obj(obj):
        if isinstance(obj, dict):
            if obj.get("sequence"):
                return obj.get("sequence", "")
            if obj.get("seq"):
                return obj.get("seq", "")
            if obj.get("primer"):
                return str(obj.get("primer", ""))
        elif isinstance(obj, str):
            return obj
        return ""

    @staticmethod
    def _iter_candidate_collections_from_result(result: dict):
        if not isinstance(result, dict):
            return []
        collections = []
        for key in ("pairs", "blast_results", "evaluated_blast_results"):
            value = result.get(key, [])
            if isinstance(value, list):
                collections.append((key, value))
        return collections

    def _get_item_by_pair_index_from_result(self, result: dict, collection_name: str, pair_index):
        if not isinstance(result, dict):
            return {}

        value = result.get(collection_name, [])
        if not isinstance(value, list):
            return {}

        target = self._normalize_pair_index(pair_index)

        for i, item in enumerate(value):
            if not isinstance(item, dict):
                continue

            item_idx = self._normalize_pair_index(item.get("pair_index"))
            if item_idx == target:
                return item
            if i == target:
                return item

        return {}

    def _find_best_pair_source_from_result(self, result: dict, pair_index):
        target = self._normalize_pair_index(pair_index)

        for key, collection in self._iter_candidate_collections_from_result(result):
            for i, item in enumerate(collection):
                if not isinstance(item, dict):
                    continue
                item_idx = self._normalize_pair_index(item.get("pair_index"))
                if item_idx == target or i == target:
                    return key, item

        return "", {}

    def _lookup_forward_from_result(self, result: dict, item: dict) -> str:
        seq = self._extract_sequence_from_obj(item.get("forward"))
        if seq:
            return seq

        for key in ("pair_profile", "pair"):
            block = item.get(key, {}) or {}
            seq = self._extract_sequence_from_obj(block.get("forward"))
            if seq:
                return seq

        pair_index = item.get("pair_index")
        for collection_name in ("pairs", "blast_results", "evaluated_blast_results"):
            source_obj = self._get_item_by_pair_index_from_result(result, collection_name, pair_index)
            if not source_obj:
                continue

            seq = self._extract_sequence_from_obj(source_obj.get("forward"))
            if seq:
                return seq

            for key in ("pair_profile", "pair"):
                block = source_obj.get(key, {}) or {}
                seq = self._extract_sequence_from_obj(block.get("forward"))
                if seq:
                    return seq

        return ""

    def _lookup_reverse_from_result(self, result: dict, item: dict) -> str:
        seq = self._extract_sequence_from_obj(item.get("reverse"))
        if seq:
            return seq

        for key in ("pair_profile", "pair"):
            block = item.get(key, {}) or {}
            seq = self._extract_sequence_from_obj(block.get("reverse"))
            if seq:
                return seq

        pair_index = item.get("pair_index")
        for collection_name in ("pairs", "blast_results", "evaluated_blast_results"):
            source_obj = self._get_item_by_pair_index_from_result(result, collection_name, pair_index)
            if not source_obj:
                continue

            seq = self._extract_sequence_from_obj(source_obj.get("reverse"))
            if seq:
                return seq

            for key in ("pair_profile", "pair"):
                block = source_obj.get(key, {}) or {}
                seq = self._extract_sequence_from_obj(block.get("reverse"))
                if seq:
                    return seq

        return ""

    def _lookup_product_size_from_result(self, result: dict, item: dict):
        v = item.get("product_size")
        if v not in (None, ""):
            return v

        pair_index = item.get("pair_index")
        for collection_name in ("pairs", "blast_results", "evaluated_blast_results"):
            source_obj = self._get_item_by_pair_index_from_result(result, collection_name, pair_index)
            if not source_obj:
                continue
            v = source_obj.get("product_size")
            if v not in (None, ""):
                return v

        return ""

    @staticmethod
    def _short_text(text: str, max_len: int = 18) -> str:
        if not text:
            return ""
        text = str(text)
        if len(text) <= max_len:
            return text
        return text[:max_len] + "..."

    # -------------------------
    # Result table
    # -------------------------
    def _fill_result_table(self, rows: list):
        self._clear_results()

        if not rows:
            self.result_hint.config(text="Run completed, but no combined ranking rows were found.")
            if self.result_data:
                self.export_json_button.config(state="normal")
            return

        for row in rows:
            values = (
                row.get("rank", ""),
                row.get("pair_index", ""),
                row.get("forward", ""),
                row.get("reverse", ""),
                row.get("product_size", ""),
                row.get("source", ""),
                row.get("score", ""),
                row.get("label", ""),
                row.get("warnings", ""),
            )
            self.tree.insert("", tk.END, values=values)

        self.result_hint.config(text=f"{len(rows)} rows shown from combined ranking.")
        self.export_json_button.config(state="normal")
        self.export_csv_button.config(state="normal")

    def _ranking_to_rows(self, result: dict, top_n: int):
        rankings = result.get("rankings", {}) or {}
        combined = rankings.get("combined_ranking", []) or []

        rows = []
        for i, item in enumerate(combined[:top_n], start=1):
            if not isinstance(item, dict):
                continue

            forward_full = self._lookup_forward_from_result(result, item)
            reverse_full = self._lookup_reverse_from_result(result, item)
            product_size = self._lookup_product_size_from_result(result, item)

            row = {
                "rank": i,
                "pair_index": item.get("pair_index", ""),
                "forward": self._short_text(forward_full, 18),
                "reverse": self._short_text(reverse_full, 18),
                "forward_full": forward_full,
                "reverse_full": reverse_full,
                "product_size": product_size,
                "source": item.get("source", item.get("tissue", item.get("site", ""))),
                "score": self._format_score(
                    item.get("combined_score", item.get("score", ""))
                ),
                "label": item.get(
                    "recommendation_label",
                    item.get("final_recommendation", item.get("label", "")),
                ),
                "warnings": self._flatten_warnings(item),
                "_raw": item,
            }
            rows.append(row)

        return rows

    @staticmethod
    def _format_score(value):
        if value in (None, ""):
            return ""
        try:
            return f"{float(value):.4f}"
        except Exception:
            return str(value)

    @staticmethod
    def _flatten_warnings(item: dict) -> str:
        parts = []
        for key in ("warnings", "source_warnings", "tags", "source_tags", "all_tags"):
            value = item.get(key)
            if isinstance(value, list):
                parts.extend(str(x) for x in value if x not in (None, ""))
            elif isinstance(value, str) and value.strip():
                parts.append(value.strip())

        seen = set()
        deduped = []
        for p in parts:
            if p not in seen:
                seen.add(p)
                deduped.append(p)

        return "; ".join(deduped)

    # -------------------------
    # Clipboard helpers
    # -------------------------
    def _copy_to_clipboard(self, text: str, label: str = "Text"):
        try:
            self.root.clipboard_clear()
            self.root.clipboard_append(text or "")
            self.root.update()
            self._append_log(f"[INFO] {label} copied to clipboard.")
        except Exception as exc:
            messagebox.showerror(APP_TITLE, f"Failed to copy {label}:\n{exc}")

    def _make_readonly_copyable_entry(self, parent, value: str, row: int, col: int, width: int = 42):
        entry = ttk.Entry(parent, width=width)
        entry.grid(row=row, column=col, sticky="ew", pady=3)
        entry.insert(0, value or "")
        entry.configure(state="readonly")
        return entry

    # -------------------------
    # Detail window
    # -------------------------
    def _on_row_double_click(self, event=None):
        selected = self.tree.selection()
        if not selected:
            return

        item_id = selected[0]
        index = self.tree.index(item_id)
        if index < 0 or index >= len(self.current_rows):
            return

        row = self.current_rows[index]
        raw = row.get("_raw", {}) or {}
        self._open_detail_window(row, raw)

    def _open_detail_window(self, row: dict, raw: dict):
        win = tk.Toplevel(self.root)
        win.title(f"Pair Detail - rank {row.get('rank', '')}")
        win.geometry("1080x780")
        win.minsize(920, 680)

        container = ttk.Frame(win, padding=12)
        container.pack(fill="both", expand=True)

        pair_index = row.get("pair_index", "")
        source_name, source_obj = self._find_best_pair_source_from_result(self.result_data, pair_index)

        forward_seq = row.get("forward_full", "") or self._extract_sequence_from_obj(
            source_obj.get("forward") if isinstance(source_obj, dict) else {}
        )
        reverse_seq = row.get("reverse_full", "") or self._extract_sequence_from_obj(
            source_obj.get("reverse") if isinstance(source_obj, dict) else {}
        )
        product_size = row.get("product_size", "") or (
            source_obj.get("product_size", "") if isinstance(source_obj, dict) else ""
        )

        summary = ttk.LabelFrame(container, text="Summary", padding=10)
        summary.pack(fill="x", expand=False, pady=(0, 10))
        for c in range(5):
            summary.columnconfigure(c, weight=1)

        ttk.Label(summary, text="Rank:", font=("Arial", 10, "bold")).grid(
            row=0, column=0, sticky="w", padx=(0, 6), pady=4
        )
        ttk.Label(summary, text=str(row.get("rank", ""))).grid(row=0, column=1, sticky="w", pady=4)

        ttk.Label(summary, text="Pair Index:", font=("Arial", 10, "bold")).grid(
            row=0, column=2, sticky="w", padx=(12, 6), pady=4
        )
        ttk.Label(summary, text=str(pair_index)).grid(row=0, column=3, sticky="w", pady=4)

        ttk.Label(summary, text="Forward:", font=("Arial", 10, "bold")).grid(
            row=1, column=0, sticky="nw", padx=(0, 6), pady=4
        )
        self._make_readonly_copyable_entry(summary, forward_seq, 1, 1, width=48)
        ttk.Button(
            summary,
            text="Copy Forward",
            command=lambda: self._copy_to_clipboard(forward_seq, "Forward primer"),
        ).grid(row=1, column=2, sticky="w", padx=(12, 0), pady=4)

        ttk.Label(summary, text="Reverse:", font=("Arial", 10, "bold")).grid(
            row=2, column=0, sticky="nw", padx=(0, 6), pady=4
        )
        self._make_readonly_copyable_entry(summary, reverse_seq, 2, 1, width=48)
        ttk.Button(
            summary,
            text="Copy Reverse",
            command=lambda: self._copy_to_clipboard(reverse_seq, "Reverse primer"),
        ).grid(row=2, column=2, sticky="w", padx=(12, 0), pady=4)

        ttk.Label(summary, text="Product Size:", font=("Arial", 10, "bold")).grid(
            row=3, column=0, sticky="w", padx=(0, 6), pady=4
        )
        ttk.Label(summary, text=str(product_size)).grid(row=3, column=1, sticky="w", pady=4)

        ttk.Label(summary, text="Source:", font=("Arial", 10, "bold")).grid(
            row=3, column=2, sticky="w", padx=(12, 6), pady=4
        )
        ttk.Label(summary, text=str(row.get("source", ""))).grid(row=3, column=3, sticky="w", pady=4)

        ttk.Label(summary, text="Score:", font=("Arial", 10, "bold")).grid(
            row=4, column=0, sticky="w", padx=(0, 6), pady=4
        )
        ttk.Label(summary, text=str(row.get("score", ""))).grid(row=4, column=1, sticky="w", pady=4)

        ttk.Label(summary, text="Label:", font=("Arial", 10, "bold")).grid(
            row=4, column=2, sticky="w", padx=(12, 6), pady=4
        )
        ttk.Label(summary, text=str(row.get("label", ""))).grid(row=4, column=3, sticky="w", pady=4)

        warn_frame = ttk.LabelFrame(container, text="Warnings", padding=10)
        warn_frame.pack(fill="x", expand=False, pady=(0, 10))
        ttk.Label(
            warn_frame,
            text=row.get("warnings", "") or "(none)",
            justify="left",
            wraplength=980,
        ).pack(anchor="w")

        source_frame = ttk.LabelFrame(container, text="Resolved Pair Source", padding=10)
        source_frame.pack(fill="x", expand=False, pady=(0, 10))
        ttk.Label(
            source_frame,
            text=f"Loaded from: {source_name or '(not found)'}",
            font=("Arial", 10, "bold"),
        ).pack(anchor="w")

        if isinstance(source_obj, dict):
            extra_lines = []
            if source_obj.get("forward_warning"):
                extra_lines.append(f"Forward warning: {source_obj.get('forward_warning')}")
            if source_obj.get("reverse_warning"):
                extra_lines.append(f"Reverse warning: {source_obj.get('reverse_warning')}")
            if source_obj.get("tm_diff") is not None:
                extra_lines.append(f"Tm diff: {source_obj.get('tm_diff')}")
            if source_obj.get("heterodimer_dG") is not None:
                extra_lines.append(f"Heterodimer dG: {source_obj.get('heterodimer_dG')}")
            if extra_lines:
                ttk.Label(
                    source_frame,
                    text="\n".join(extra_lines),
                    justify="left",
                    wraplength=980,
                ).pack(anchor="w", pady=(6, 0))

        raw_frame = ttk.LabelFrame(container, text="Detail JSON", padding=10)
        raw_frame.pack(fill="both", expand=True)

        tools = ttk.Frame(raw_frame)
        tools.pack(fill="x", pady=(0, 8))

        detail_payload = {
            "ranking_item": raw,
            "resolved_pair_source_name": source_name,
            "resolved_pair_source_item": source_obj,
        }
        try:
            pretty = json.dumps(detail_payload, ensure_ascii=False, indent=2)
        except Exception:
            pretty = str(detail_payload)

        ttk.Button(
            tools,
            text="Copy All JSON",
            command=lambda: self._copy_to_clipboard(pretty, "Detail JSON"),
        ).pack(side="right")

        text = ScrolledText(raw_frame, wrap=tk.WORD, font=("Consolas", 10))
        text.pack(fill="both", expand=True)
        text.insert("1.0", pretty)
        text.configure(state="normal")
        text.bind("<Key>", lambda e: "break")

    # -------------------------
    # Export
    # -------------------------
    def _on_export_json(self):
        if not self.result_data:
            messagebox.showinfo(APP_TITLE, "No result data to export.")
            return

        path = filedialog.asksaveasfilename(
            title="Export JSON",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(self.result_data, f, ensure_ascii=False, indent=2)
            self._append_log(f"[INFO] JSON exported to: {path}")
        except Exception as exc:
            messagebox.showerror(APP_TITLE, f"Failed to export JSON:\n{exc}")

    def _on_export_csv(self):
        if not self.current_rows:
            messagebox.showinfo(APP_TITLE, "No table rows to export.")
            return

        path = filedialog.asksaveasfilename(
            title="Export CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            fieldnames = [
                "rank",
                "pair_index",
                "forward",
                "reverse",
                "product_size",
                "source",
                "score",
                "label",
                "warnings",
            ]
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for row in self.current_rows:
                    writer.writerow({k: row.get(k, "") for k in fieldnames})
            self._append_log(f"[INFO] CSV exported to: {path}")
        except Exception as exc:
            messagebox.showerror(APP_TITLE, f"Failed to export CSV:\n{exc}")


def main():
    root = tk.Tk()
    try:
        style = ttk.Style()
        if "vista" in style.theme_names():
            style.theme_use("vista")
        elif "clam" in style.theme_names():
            style.theme_use("clam")
    except Exception:
        pass

    PrimerScoreGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()