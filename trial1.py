# trial1.py - COMPLETE VERSION WITH INTEGRATED SOLVER
import tkinter as tk
from tkinter import ttk, simpledialog, messagebox
from PIL import Image, ImageTk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import math
import copy
import sys
import ctypes
import os
import re
import networkx as nx
import numpy as np

# IMPORT THE NUMERICAL SOLVER
import graphz

# --- Existing Graph Logic ---


def get_base_path():
    """Returns the base path for assets, handling PyInstaller's temporary folder."""
    # sys._MEIPASS is the path created by PyInstaller for temporary files
    if getattr(sys, 'frozen', False):
        return sys._MEIPASS

    # Path when running directly from the source folder
    return os.path.dirname(os.path.abspath(__file__))


def generate_graph_from_netlist(filename="output.txt"):
    G = nx.MultiGraph()
    if not os.path.exists(filename):
        return None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r'\s+', line)
            if len(parts) >= 3:
                name, node1, node2 = parts[0], parts[1], parts[2]
                G.add_edge(node1, node2, label=name)
    return G

# --- Main GUI Class ---


class CircuitGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("ARHSpice_Solver")
        self.base_path = get_base_path()

        if sys.platform == 'win32':
            try:
                ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
                    u'mycompany.arhspice.tabs')
            except:
                pass

        self.geometry("1200x800")

        # --- STORAGE ---
        self.component_list = []
        self.volt_vars = {}
        self.curr_vars = {}
        self.time_entries = {}

        # --- TAB SYSTEM ---
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.editor_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.editor_frame, text="  Schematic Editor  ")

        self.voltage_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.voltage_frame, text="  Voltage Graphs  ")
        self.setup_result_tab(self.voltage_frame, "voltage")

        self.current_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.current_frame, text="  Current Graphs  ")
        self.setup_result_tab(self.current_frame, "current")

        # --- EDITOR UI ---
        self.setup_editor_ui()

        # Logic Init
        self.circuit = CircuitGraph()
        self.component_counters = {
            "R": 1, "C": 1, "L": 1, "V": 1, "I": 1, "G": 1}
        self.selected_component_type, self.wire_mode, self.delete_mode = None, False, False
        self.selected_terminals, self._drag_data = [], {
            "x": 0, "y": 0, "item": None}
        self.unit_map, self.history_stack = {
            'R': 'Ω', 'C': 'F', 'L': 'H', 'V': 'V', 'I': 'A'}, []

        # Bindings
        self.canvas.bind("<Button-1>", self.on_canvas_click)
        self.canvas.tag_bind("component", "<ButtonPress-1>",
                             self.on_component_press)
        self.canvas.tag_bind("component", "<B1-Motion>",
                             self.on_component_drag)
        self.canvas.tag_bind(
            "component", "<ButtonRelease-1>", self.on_component_release)
        self.canvas.tag_bind("component", "<Button-3>", self.rotate_component)
        self.canvas.tag_bind(
            "component", "<Double-Button-1>", self.edit_component_value)
        self.bind_all("<Control-z>", self.undo_action)
        self.update_status_text("Component Placement Mode")
        self.save_state()

    def setup_editor_ui(self):
        style = ttk.Style(self)
        style.configure("TButton", padding=5, font=('Arial', 9))
        style.configure("Control.TFrame", background='#ECECEC')
        style.configure("Green.TButton", background="#D5E8D4",
                        foreground="black")
        style.configure("Red.TButton", background="#F8CECC",
                        foreground="black")
        style.configure("Blue.TButton", background="#DAE8FC",
                        foreground="black")

        control_frame = ttk.Frame(
            self.editor_frame, style="Control.TFrame", padding=(10, 5))
        control_frame.pack(side=tk.TOP, fill=tk.X)

        self.canvas = tk.Canvas(
            self.editor_frame, bg="white", highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)

        self.component_images, self.photo_images = {}, {}
        if not self.load_images():
            self.destroy()
            return

        comp_buttons = {"Resistor": self.select_resistor, "Capacitor": self.select_capacitor,
                        "Inductor": self.select_inductor, "Voltage Source": self.select_voltage_source,
                        "Current Source": self.select_current_source, "Ground": self.select_ground}

        for text, command in comp_buttons.items():
            ttk.Button(control_frame, text=text, command=command,
                       style="TButton").pack(side=tk.LEFT, padx=3)

        ttk.Separator(control_frame, orient='vertical').pack(
            side=tk.LEFT, padx=10, fill='y')
        ttk.Button(control_frame, text="Wire",
                   command=self.toggle_wire_mode).pack(side=tk.LEFT, padx=3)
        ttk.Button(control_frame, text="Delete", command=self.toggle_delete_mode,
                   style="Red.TButton").pack(side=tk.LEFT, padx=3)
        ttk.Button(control_frame, text="Undo", command=self.undo_action,
                   style="Blue.TButton").pack(side=tk.LEFT, padx=3)

        ttk.Button(control_frame, text="▶ SIMULATE", command=self.simulate,
                   style="Green.TButton").pack(side=tk.RIGHT, padx=5)

    def setup_result_tab(self, parent_frame, plot_type):
        container = ttk.Frame(parent_frame)
        container.pack(fill=tk.BOTH, expand=True)

        sidebar = ttk.Frame(container, width=250, padding=10)
        sidebar.pack(side=tk.LEFT, fill=tk.Y)

        ttk.Label(sidebar, text=f"{plot_type.title()} Settings", font=(
            "Arial", 12, "bold")).pack(pady=5)

        time_frame = ttk.LabelFrame(
            sidebar, text="Time Interval (seconds)", padding=5)
        time_frame.pack(fill=tk.X, pady=10)

        tk.Label(time_frame, text="Start (t1):").grid(row=0, column=0, padx=2)
        t1_entry = ttk.Entry(time_frame, width=8)
        t1_entry.insert(0, "0")
        t1_entry.grid(row=0, column=1, padx=2)

        tk.Label(time_frame, text="End (t2):").grid(
            row=1, column=0, padx=2, pady=5)
        t2_entry = ttk.Entry(time_frame, width=8)
        t2_entry.insert(0, "0.1")
        t2_entry.grid(row=1, column=1, padx=2, pady=5)

        self.time_entries[plot_type] = (t1_entry, t2_entry)

        # BUTTON
        ttk.Button(sidebar, text="Update Plot", style="Green.TButton",
                   command=lambda: self.update_plot(plot_type)).pack(fill=tk.X, pady=10)

        comp_frame = ttk.LabelFrame(
            sidebar, text="Select Components", padding=5)
        comp_frame.pack(fill=tk.BOTH, expand=True, pady=10)

        canvas = tk.Canvas(comp_frame, height=300)
        scrollbar = ttk.Scrollbar(
            comp_frame, orient="vertical", command=canvas.yview)
        scroll_frame = ttk.Frame(canvas)

        scroll_frame.bind("<Configure>", lambda e: canvas.configure(
            scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        if plot_type == "voltage":
            self.v_chk_frame = scroll_frame
        else:
            self.c_chk_frame = scroll_frame

        graph_area = ttk.Frame(container, style="Control.TFrame")
        graph_area.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        ttk.Label(graph_area, text="Select components and click 'Update Plot'").place(
            relx=0.5, rely=0.5, anchor="center")

        if plot_type == "voltage":
            self.v_graph_frame = graph_area
        else:
            self.c_graph_frame = graph_area

    # --- Loading & Helper Methods ---
    def load_images(self):
        component_files = {"resistor": "resistor.png", "capacitor": "capacitor.png", "inductor": "inductor.png",
                           "voltage_source": "voltage_source.png", "current_source": "current_source.png",
                           "ground": "ground.png"}
        try:
            for name, filename in component_files.items():
                full_path = os.path.join(self.base_path, filename)
                img = Image.open(full_path)
                if name == "ground":
                    img = img.convert("RGBA")
                if name in ["voltage_source", "current_source"]:
                    size = (40, 60)
                elif name == "ground":
                    size = (45, 45)
                else:
                    size = (120, 30)
                resample = Image.Resampling.LANCZOS if hasattr(
                    Image.Resampling, 'LANCZOS') else Image.ANTIALIAS
                self.component_images[name] = img.resize(size, resample)
            return True
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load images: {e}")
            return False

    def save_state(self):
        state = copy.deepcopy(
            {'components': self.circuit.components, 'connections': self.circuit.connections})
        self.history_stack.append(state)
        if len(self.history_stack) > 30:
            self.history_stack.pop(0)

    def undo_action(self, event=None):
        if len(self.history_stack) > 1:
            self.history_stack.pop()
            prev = self.history_stack[-1]
            self.circuit.components = copy.deepcopy(prev['components'])
            self.circuit.connections = copy.deepcopy(prev['connections'])
            self.redraw_canvas()

    def redraw_canvas(self):
        self.canvas.delete("all")
        self.photo_images.clear()
        for cid, comp in self.circuit.components.items():
            ctype = self.get_type_from_id(cid)
            img = self.component_images[ctype].rotate(-comp.angle, expand=True)
            self.photo_images[cid] = ImageTk.PhotoImage(img)
            self.canvas.create_image(
                comp.x, comp.y, image=self.photo_images[cid], tags=(cid, "component"))
            self.draw_terminals(comp)
            self.update_component_text(cid)
        for conn in self.circuit.connections:
            c1, c2 = self.circuit.components.get(
                conn['c1']), self.circuit.components.get(conn['c2'])
            if c1 and c2:
                x1, y1 = c1.terminals[conn['t1']]
                x2, y2 = c2.terminals[conn['t2']]
                conn['wire_id'] = self.canvas.create_line(
                    x1, y1, x2, y2, fill="black", width=2)
                self.canvas.tag_lower(conn['wire_id'], "component")

    # --- Mode Setters ---
    def set_mode(self, comp_type=None, wire=False, delete=False):
        self.selected_component_type, self.wire_mode, self.delete_mode = comp_type, wire, delete
        self.canvas.config(cursor="X_cursor" if delete else "")
        self.update_status_text(
            f"Mode: {comp_type if comp_type else ('Wire' if wire else ('Delete' if delete else 'Standard'))}")

    def select_resistor(self): self.set_mode(comp_type="resistor")
    def select_capacitor(self): self.set_mode(comp_type="capacitor")
    def select_inductor(self): self.set_mode(comp_type="inductor")
    def select_voltage_source(self): self.set_mode(comp_type="voltage_source")
    def select_current_source(self): self.set_mode(comp_type="current_source")
    def select_ground(self): self.set_mode(comp_type="ground")
    def toggle_wire_mode(self): self.set_mode(wire=not self.wire_mode)
    def toggle_delete_mode(self): self.set_mode(delete=not self.delete_mode)

    def update_status_text(self, text):
        self.canvas.delete("status_text")
        self.canvas.create_text(10, 10, text=text, tags="status_text", anchor="nw", font=(
            "Arial", 9), fill="gray30")

    # --- Interactions ---
    def on_canvas_click(self, event):
        if self.selected_component_type and not (self.wire_mode or self.delete_mode):
            self.place_component(event.x, event.y)

    def place_component(self, x, y):
        self.save_state()
        ctype = self.selected_component_type
        # Ensure 'I' is used for Current Sources
        prefix = "L" if ctype == "inductor" else (
            "I" if ctype == "current_source" else ctype[0].upper())
        cid = f"{prefix}{self.component_counters[prefix]}"
        self.component_counters[prefix] += 1

        val = None
        if ctype != "ground":
            val = self.validate_and_get_value(ctype, f"Value for {cid}:")
            if val is None:
                self.component_counters[prefix] -= 1
                self.undo_action()
                return

        self.photo_images[cid] = ImageTk.PhotoImage(
            self.component_images[ctype])
        self.canvas.create_image(
            x, y, image=self.photo_images[cid], tags=(cid, "component"))
        terminals = self.calculate_terminal_positions(cid, x, y, 0)
        comp = Component(cid, terminals, val, x, y, 0)
        self.circuit.add_component(comp)
        self.draw_terminals(comp)
        self.update_component_text(cid)
        self.set_mode()

    def validate_and_get_value(self, ctype, prompt, init=""):
        while True:
            val = simpledialog.askstring(
                f"{ctype} Value", prompt, initialvalue=init)
            if val is None:
                return None
            if val:
                return val
            messagebox.showerror("Error", "Invalid Input")

    def edit_component_value(self, event):
        try:
            cid = self.canvas.gettags(
                self.canvas.find_closest(event.x, event.y)[0])[0]
        except:
            return
        comp = self.circuit.components.get(cid)
        if not comp or cid.startswith("G"):
            return
        self.save_state()
        new_val = self.validate_and_get_value(
            self.get_type_from_id(cid), "New Value:", comp.value)
        if new_val:
            comp.value = new_val
            self.update_component_text(cid)
        else:
            self.undo_action()

    def update_component_text(self, cid):
        self.canvas.delete(f"text_{cid}")
        comp = self.circuit.components.get(cid)
        if comp and comp.value:
            unit = self.unit_map.get(cid[0], "")
            txt = f"{comp.value}{unit}"
            angle = comp.angle % 360
            xo, yo = 0, 0
            if angle == 0:
                yo = 22
            elif angle == 90:
                xo = 40
            elif angle == 180:
                yo = -22
            elif angle == 270:
                xo = -40
            self.canvas.create_text(comp.x + xo, comp.y + yo, text=txt, tags=(
                f"text_{cid}", "component_text"), fill="navy", font=("Arial", 8, "bold"))

    def calculate_terminal_positions(self, cid, x, y, angle):
        ctype = self.get_type_from_id(cid)
        if cid.startswith("G"):
            return {"Ground": (x, y - 20)}
        w, _ = self.component_images[ctype].size
        offset = (w/2) * 0.92
        rad = math.radians(angle)
        return {f"{cid}.n1": (x - offset*math.cos(rad), y - offset*math.sin(rad)),
                f"{cid}.n2": (x + offset*math.cos(rad), y + offset*math.sin(rad))}

    def draw_terminals(self, comp):
        self.canvas.delete(f"terminals_{comp.component_id}")
        for t, (tx, ty) in comp.terminals.items():
            tag = f"term_{t}"
            self.canvas.create_oval(
                tx-4, ty-4, tx+4, ty+4, fill="black", tags=(tag, f"terminals_{comp.component_id}"))
            self.canvas.tag_bind(tag, "<Button-1>", lambda e,
                                 c=comp.component_id, t=t: self.select_terminal(c, t))

    def select_terminal(self, cid, term):
        if not self.wire_mode:
            return
        if not self.selected_terminals:
            self.selected_terminals.append((cid, term))
        else:
            (c1, t1) = self.selected_terminals[0]
            if c1 == cid:
                return
            self.save_state()
            x1, y1 = self.circuit.components[c1].terminals[t1]
            x2, y2 = self.circuit.components[cid].terminals[term]
            wid = self.canvas.create_line(
                x1, y1, x2, y2, fill="black", width=2)
            self.canvas.tag_lower(wid, "component")
            self.circuit.add_connection(c1, t1, cid, term, wid)
            self.selected_terminals = []

    def on_component_press(self, event):
        item = self.canvas.find_closest(event.x, event.y)[0]
        cid = self.canvas.gettags(item)[0]
        if self.delete_mode:
            self.delete_component(cid)
            return
        if not self.wire_mode:
            self.save_state()
            self._drag_data.update({"item": item, "x": event.x, "y": event.y})

    def on_component_drag(self, event):
        if self.wire_mode or self.delete_mode or not self._drag_data["item"]:
            return
        dx, dy = event.x - self._drag_data["x"], event.y - self._drag_data["y"]
        cid = self.canvas.gettags(self._drag_data["item"])[0]
        self.canvas.move(self._drag_data["item"], dx, dy)
        self.canvas.move(f"text_{cid}", dx, dy)
        self._drag_data.update({"x": event.x, "y": event.y})
        comp = self.circuit.components[cid]
        comp.x += dx
        comp.y += dy
        comp.terminals = self.calculate_terminal_positions(
            cid, comp.x, comp.y, comp.angle)
        self.draw_terminals(comp)
        self.update_connected_wires(cid)

    def on_component_release(self, event): self._drag_data["item"] = None

    def rotate_component(self, event):
        self.save_state()
        try:
            cid = self.canvas.gettags(
                self.canvas.find_closest(event.x, event.y)[0])[0]
        except:
            return
        comp = self.circuit.components[cid]
        comp.angle = (comp.angle + 90) % 360
        img = self.component_images[self.get_type_from_id(
            cid)].rotate(-comp.angle, expand=True)
        self.photo_images[cid] = ImageTk.PhotoImage(img)
        self.canvas.itemconfig(cid, image=self.photo_images[cid])
        comp.terminals = self.calculate_terminal_positions(
            cid, comp.x, comp.y, comp.angle)
        self.draw_terminals(comp)
        self.update_connected_wires(cid)
        self.update_component_text(cid)

    def delete_component(self, cid):
        self.save_state()
        self.canvas.delete(cid)
        self.canvas.delete(f"terminals_{cid}")
        self.canvas.delete(f"text_{cid}")
        for conn in [c for c in self.circuit.connections if c['c1'] == cid or c['c2'] == cid]:
            self.canvas.delete(conn['wire_id'])
            self.circuit.connections.remove(conn)
        if cid in self.circuit.components:
            del self.circuit.components[cid]

    def update_connected_wires(self, cid):
        for conn in [c for c in self.circuit.connections if c['c1'] == cid or c['c2'] == cid]:
            c1, c2 = self.circuit.components.get(
                conn['c1']), self.circuit.components.get(conn['c2'])
            if c1 and c2:
                x1, y1 = c1.terminals[conn['t1']]
                x2, y2 = c2.terminals[conn['t2']]
                self.canvas.coords(conn['wire_id'], x1, y1, x2, y2)

    def get_type_from_id(self, cid):
        return {'R': 'resistor', 'C': 'capacitor', 'L': 'inductor', 'V': 'voltage_source', 'I': 'current_source', 'G': 'ground'}.get(cid[0])

    # --- SIMULATION LOGIC ---
    def simulate(self):
        if not self.circuit.components:
            messagebox.showwarning("Empty", "No components found.")
            return

        # --- GROUND CHECK ---
        has_ground = False
        for cid, comp in self.circuit.components.items():
            if cid.startswith("G"):
                has_ground = True
                break
        if not has_ground:
            messagebox.showerror(
                "Circuit Error", "Missing GROUND.\nPlease place a Ground component.")
            return

        netlist = self.circuit.generate_netlist()
        netlist_str = ""
        try:
            with open("output.txt", "w") as file:
                for entry in netlist:
                    line = " ".join(map(str, entry))
                    file.write(line + "\n")
                    netlist_str += line + "\n"
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save netlist: {e}")
            return

        # 1. PARSE STRUCTURE
        self.component_list = graphz.parse_netlist_structure(netlist_str)
        self.update_ui_options()

        # 2. Switch tab
        self.notebook.select(self.voltage_frame)

        # 3. Show success message
        messagebox.showinfo("Simulation Ready",
                            "Circuit netlist generated successfully!\n"
                            "Switch to Voltage/Current tabs to view results.")

    def update_ui_options(self):
        """Creates checkboxes with values"""
        def populate_frame(frame, var_dict):
            for widget in frame.winfo_children():
                widget.destroy()
            var_dict.clear()

            for b in self.component_list:
                name = b['name']
                val_str = b.get('val_str', '')
                text_label = f"{name} ({val_str})"

                var = tk.BooleanVar(value=False)
                var_dict[name] = var
                cb = ttk.Checkbutton(frame, text=text_label, variable=var)
                cb.pack(anchor="w", padx=5)

        populate_frame(self.v_chk_frame, self.volt_vars)
        populate_frame(self.c_chk_frame, self.curr_vars)

    def update_plot(self, plot_type):
        """Triggers the actual simulation logic in graphz"""
        if not self.component_list:
            messagebox.showinfo("Info", "Please run 'Simulate' first.")
            return

        # 1. Get Time Intervals
        t_start_entry, t_end_entry = self.time_entries[plot_type]
        try:
            t1 = float(t_start_entry.get())
            t2 = float(t_end_entry.get())
            if t1 < 0 or t2 <= t1:
                raise ValueError
        except:
            messagebox.showerror(
                "Input Error", "Invalid time interval.\nEnsure End Time > Start Time >= 0")
            return

        # 2. Get Selected Components
        var_dict = self.volt_vars if plot_type == "voltage" else self.curr_vars
        selected_names = [name for name, var in var_dict.items() if var.get()]

        if not selected_names:
            messagebox.showinfo(
                "Selection", "Please select at least one component to plot.")
            return

        # 3. RUN SIMULATION using the integrated Python solver
        try:
            with open("output.txt", "r") as f:
                netlist_str = f.read()

            # CALL THE INTEGRATED PYTHON SOLVER (replaces MATLAB)
            results = graphz.solve_circuit_mna(netlist_str, t1, t2)

            if isinstance(results, str):
                # Error message returned string
                messagebox.showerror("Simulation Error", results)
                return

            if results is None:
                messagebox.showerror("Simulation Error",
                                     "Simulation failed. Check connectivity.")
                return

            # 4. Filter and Plot
            filtered_data = [
                res for res in results if res['name'] in selected_names]

            frame = self.v_graph_frame if plot_type == "voltage" else self.c_graph_frame
            color = 'orange' if plot_type == "voltage" else 'blue'

            self.generate_custom_plot(frame, filtered_data, plot_type, color)

        except Exception as e:
            messagebox.showerror(
                "Solver Error", f"Error during calculation: {e}")

    def generate_custom_plot(self, frame, results_data, plot_key, color):
        for widget in frame.winfo_children():
            widget.destroy()

        num = len(results_data)
        if num == 0:
            return

        fig = Figure(figsize=(5, 4), dpi=100)

        if num == 1:
            rows, cols = 1, 1
        elif num == 2:
            rows, cols = 2, 1
        elif num <= 4:
            rows, cols = 2, 2
        else:
            rows, cols = math.ceil(num/2), 2

        for i, res in enumerate(results_data):
            ax = fig.add_subplot(rows, cols, i+1)
            name = res['name']
            val = res['val_str']
            time = res['time']
            # Select Voltage or Current
            y_vals = res['voltage'] if plot_key == 'voltage' else res['current']

            ax.plot(time, y_vals, color=color)
            ax.set_title(f"{name} ({val})")
            ax.grid(True, alpha=0.5)

        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, frame)
        toolbar.update()

        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# --- Helper Classes ---


class DisjointSet:
    def __init__(self): self.parent = {}; self.rank = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        rootX, rootY = self.find(x), self.find(y)
        if rootX != rootY:
            if rootY == 'Ground':
                self.parent[rootX] = rootY
                return
            if rootX == 'Ground':
                self.parent[rootY] = rootX
                return
            self.parent[rootY] = rootX

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x


class Component:
    def __init__(self, component_id, terminals, value=None, x=0, y=0, angle=0):
        self.component_id, self.terminals, self.value = component_id, terminals, value
        self.x, self.y, self.angle = x, y, angle


class CircuitGraph:
    def __init__(self): self.components = {}; self.connections = []

    def add_component(
        self, component): self.components[component.component_id] = component

    def add_connection(self, c1, t1, c2, t2, wire_id):
        self.connections.append(
            {"c1": c1, "t1": t1, "c2": c2, "t2": t2, "wire_id": wire_id})

    def generate_netlist(self):
        dsu = DisjointSet()
        dsu.add('Ground')  # Abstract Ground Node

        # Register nodes
        for cid, comp in self.components.items():
            for t in comp.terminals:
                node_name = f"{cid}_{t}"
                dsu.add(node_name)
                # CRITICAL FIX: Explicitly link Ground components to the abstract 'Ground'
                if cid.startswith("G"):
                    dsu.union(node_name, 'Ground')

        # Apply wires
        for conn in self.connections:
            t1 = f"{conn['c1']}_{conn['t1']}"
            t2 = f"{conn['c2']}_{conn['t2']}"
            dsu.add(t1)
            dsu.add(t2)
            dsu.union(t1, t2)

        # Map to integers
        root_map = {}
        counter = 1

        # Find the root for ground and assign 0
        ground_root = dsu.find('Ground')
        root_map[ground_root] = 0

        # Assign numbers to others
        for item in dsu.parent:
            root = dsu.find(item)
            if root not in root_map:
                root_map[root] = counter
                counter += 1

        # Build List
        netlist = []
        for cid, comp in self.components.items():
            if cid.startswith("G"):
                continue

            # Get terminals
            terms = list(comp.terminals.keys())
            # Ensure consistent order (n1, n2)
            terms.sort()  # n1, n2

            n1_name = f"{cid}_{terms[0]}"
            n2_name = f"{cid}_{terms[1]}"

            # Get node numbers
            val1 = root_map[dsu.find(n1_name)]
            val2 = root_map[dsu.find(n2_name)]

            netlist.append([cid, val1, val2, comp.value])

        return netlist


if __name__ == "__main__":
    app = CircuitGUI()
    app.mainloop()
