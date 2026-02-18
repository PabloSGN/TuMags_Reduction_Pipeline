#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SUPER-LINTER EXTENDIDO PARA DETECTAR MODIFICACIONES NO DESEADAS EN VISTAS NUMPY
-------------------------------------------------------------------------------

Funciones nuevas:
 - Procesa múltiples archivos o carpetas recursivamente
 - Detección de aliasing entre variables
 - Análisis superficial de funciones para detectar modificaciones in-place
 - Informe HTML global con navegación por archivos
 - Colores en terminal

Autor: ChatGPT (para David Orozco)
"""

import ast
import sys
from pathlib import Path
import html

# ============================================
# Configuración: funciones modificadoras conocidas
# ============================================

KNOWN_MODIFIERS = {
    "correct_image_fft",
    "filter_frecuencies",
    "destretch",
    "balance",
    "apply_transform",
    "demodulate",
    "restore_ima",
}

# Colores terminal
RED   = "\033[91m"
YELLOW = "\033[93m"
GREEN = "\033[92m"
RESET = "\033[0m"

# ============================================
# Analizador de alias entre variables
# ============================================

class AliasTable:
    """ Gestiona alias entre variables (v1 apunta a v2, etc.) """
    def __init__(self):
        # mapa var -> names que comparten referencia
        self.alias = {}

    def add_alias(self, a, b):
        """a = b  =>  a y b son alias"""
        if a not in self.alias:
            self.alias[a] = set([a])
        if b not in self.alias:
            self.alias[b] = set([b])

        # unir conjuntos
        group = self.alias[a] | self.alias[b]
        for name in group:
            self.alias[name] = set(group)

    def is_alias(self, name):
        return name in self.alias

    def get_group(self, name):
        return self.alias.get(name, set([name]))


# ============================================
# LINTER
# ============================================

class ViewLinter(ast.NodeVisitor):

    def __init__(self, filename, global_alias):
        self.filename = filename
        self.lines = Path(filename).read_text().splitlines()
        self.views = {}          # var_name -> line number
        self.warnings = []       # tuples (line, level, msg)
        self.aliases = global_alias  # tabla compartida entre archivos

    # ---------------------------------------
    def _record(self, lineno, level, msg):
        self.warnings.append((lineno, level, msg))
        colored = msg
        if level == "ERROR":
            colored = RED + msg + RESET
        elif level == "WARN":
            colored = YELLOW + msg + RESET
        print(f"{colored}  (línea {lineno} en {self.filename})")

    # ---------------------------------------
    def mark_view(self, name, node):
        """Marca variable como vista, y extiende a alias."""
        if name is None:
            return

        # Marcar variable
        self.views[name] = node.lineno
        self._record(node.lineno, "INFO",
                     f"'{name}' es VISTA (slicing/indexación)")

        # También marcar alias existentes
        if self.aliases.is_alias(name):
            for v in self.aliases.get_group(name):
                self.views[v] = node.lineno

    # ---------------------------------------
    def extract_name(self, node):
        """Intento de extraer nombre del elemento base arr[...]"""
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Attribute):
            return node.attr
        if isinstance(node, ast.Subscript):
            return self.extract_name(node.value)
        return None

    # ========================================
    # 1) ASIGNACIONES
    # ========================================

    def visit_Assign(self, node):
        """
        Detecta:
            out = arr[...]
            out = hdul[0].data
            out = other_var     (alias)
        """
        # ---- detectar alias ----
        if len(node.targets) == 1 and isinstance(node.value, ast.Name):
            left = node.targets[0]
            if isinstance(left, ast.Name):
                self.aliases.add_alias(left.id, node.value.id)

        # ---- detectar slicing → vista ----
        if isinstance(node.value, ast.Subscript):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name):
                    self.mark_view(tgt.id, node)

        # ---- detectar FITS.data → vista ----
        if isinstance(node.value, ast.Attribute):
            if node.value.attr == "data":
                for tgt in node.targets:
                    if isinstance(tgt, ast.Name):
                        self.mark_view(tgt.id, node)

        self.generic_visit(node)

    # ========================================
    # 2) OPERACIONES IN-PLACE
    # ========================================

    def visit_AugAssign(self, node):
        """
        out += x
        out *= x
        """
        if isinstance(node.target, ast.Name):
            name = node.target.id

            # Alias: si b es alias de a, incluirlo
            group = self.aliases.get_group(name)

            for v in group:
                if v in self.views:
                    self._record(node.lineno, "ERROR",
                                 f"Modificación IN-PLACE sobre vista/alias: '{name} {node.op.__class__.__name__}='")
        self.generic_visit(node)

    # ========================================
    # 3) Escritura a vistas
    # ========================================

    def visit_Subscript(self, node):
        """
        Detecta: out[...] = ...
        """
        parent = getattr(node, "parent", None)
        if isinstance(parent, ast.Assign):
            for tgt in parent.targets:
                if tgt is node:
                    name = self.extract_name(node)
                    if name:
                        group = self.aliases.get_group(name)
                        for v in group:
                            if v in self.views:
                                self._record(node.lineno, "ERROR",
                                             f"Asignación a vista/alias: '{v}[...] = ...'")
        self.generic_visit(node)

    # ========================================
    # 4) Funciones con efectos secundarios
    # ========================================

    def visit_Call(self, node):
        func_name = None
        if isinstance(node.func, ast.Name):
            func_name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            func_name = node.func.attr

        # función potencialmente modificadora
        if func_name in KNOWN_MODIFIERS:
            for arg in node.args:
                base = self.extract_name(arg)
                if not base:
                    continue
                group = self.aliases.get_group(base)
                for v in group:
                    if v in self.views:
                        self._record(node.lineno, "WARN",
                                     f"Vista/alias '{v}' pasada a función MODIFICADORA '{func_name}'")

        self.generic_visit(node)


# ===========================
# Generador HTML Global
# ===========================

def write_global_html_report(reports, output="linter_report.html"):
    html_lines = []
    html_lines.append("<html><body>")
    html_lines.append("<h1>Informe Global — NumPy View Super-Linter</h1>")

    for filename, (lines, warnings) in reports.items():
        html_lines.append(f"<h2>{filename}</h2><pre>")
        for i, line in enumerate(lines, start=1):
            esc = html.escape(line)

            # Buscar warnings de esta línea
            w = [w for w in warnings if w[0] == i]

            if w:
                html_lines.append(
                    f'<span style="background:#ffcccc">{i:5d}: {esc}</span>'
                )
                for _, level, msg in w:
                    color = "red" if level == "ERROR" else "orange"
                    html_lines.append(f'<span style="color:{color}; font-weight:bold">    → {html.escape(msg)}</span>')
            else:
                html_lines.append(f"{i:5d}: {esc}")
        html_lines.append("</pre>")

    html_lines.append("</body></html>")

    Path(output).write_text("\n".join(html_lines))
    print(GREEN + f"\n✔️ Informe HTML creado: {output}" + RESET)


# ===========================
# Procesar archivos o carpetas
# ===========================

def collect_python_files(path):
    path = Path(path)
    if path.is_file() and path.suffix == ".py":
        return [path]
    elif path.is_dir():
        return list(path.rglob("*.py"))
    else:
        return []


# ===========================
# MAIN
# ===========================

def run(files):
    global_alias = AliasTable()
    reports = {}

    for file in files:
        print(f"\n🔎 Analizando {file}")
        src = Path(file).read_text()
        tree = ast.parse(src)

        # Añadir padres
        for n in ast.walk(tree):
            for c in ast.iter_child_nodes(n):
                c.parent = n

        # Analizar archivo
        linter = ViewLinter(file, global_alias)
        linter.visit(tree)

        reports[file] = (linter.lines, linter.warnings)

    # Informe HTML global
    write_global_html_report(reports)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python numpy_view_superlinter_v2.py archivo.py [otro.py | carpeta]")
        sys.exit(1)

    input_paths = sys.argv[1:]
    all_files = []
    for p in input_paths:
        all_files.extend(collect_python_files(p))

    run(all_files)