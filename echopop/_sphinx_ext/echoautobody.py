from sphinx.util.docutils import SphinxDirective
from sphinx.ext.autodoc import FunctionDocumenter
from docutils import nodes
import importlib
import inspect
import ast
import textwrap

# ----------------------------
# FunctionBodyDirective
# ----------------------------
class FunctionBodyDirective(SphinxDirective):
    required_arguments = 1  # function path

    def run(self):
        func_path = self.arguments[0]
        try:
            *mod_parts, func_name = func_path.split(".")
            mod_name = ".".join(mod_parts)
            module = importlib.import_module(mod_name)
            func = getattr(module, func_name)
        except Exception as e:
            return [nodes.error(text=f"Could not import function '{func_path}': {e}")]

        src_without_doc = self._get_source(func)
        return [self._make_toggle(func_name, src_without_doc)]

    def _get_source(self, func):
        try:
            src = inspect.getsource(func)
            parsed = ast.parse(src)
            func_node = parsed.body[0]

            # remove docstring
            if (len(func_node.body) > 0) and isinstance(func_node.body[0], ast.Expr) and isinstance(func_node.body[0].value, ast.Constant) and isinstance(func_node.body[0].value.value, str):
                func_node.body.pop(0)

            src_without_doc = ast.unparse(func_node)
            return textwrap.dedent(src_without_doc)
        except Exception as e:
            return f"# Could not get source: {e}"

    def _make_toggle(self, title_text, code):
        container = nodes.container(classes=["toggle", "toggle-button"])
        container['data-toggle-label-show'] = "Show source code"
        container['data-toggle-label-hide'] = "Hide source code"

        title = nodes.paragraph(text=f"Source code for {title_text}", classes=["toggle-button-title"])
        literal = nodes.literal_block(text=code, language="python", classes=["toggle-content"])

        container += title
        container += literal
        return container

class AllFunctionsDirective(SphinxDirective):
    required_arguments = 1  # full module path, e.g. echopop.module.submodule

    def run(self):
        module_path = self.arguments[0]

        try:
            module = importlib.import_module(module_path)
        except Exception as e:
            return [nodes.error(text=f"Could not import module '{module_path}': {e}")]

        output_nodes = []

        # Iterate over all functions in the module
        for name, func in inspect.getmembers(module, inspect.isfunction):
            try:
                src = inspect.getsource(func)
                parsed = ast.parse(src)
                func_node = parsed.body[0]

                # Remove docstring
                if (len(func_node.body) > 0) and isinstance(func_node.body[0], ast.Expr) and isinstance(func_node.body[0].value, ast.Constant) and isinstance(func_node.body[0].value.value, str):
                    func_node.body.pop(0)

                src_without_doc = ast.unparse(func_node)
                src_without_doc = textwrap.dedent(src_without_doc)
            except Exception as e:
                src_without_doc = f"# Could not get source for {name}: {e}"

            # Create toggle container
            container = nodes.container(classes=["toggle", "toggle-button"])
            # Optional per-container labels (may not work due to sphinx-togglebutton JS limitations)
            container['data-toggle-label-show'] = "Show source code"
            container['data-toggle-label-hide'] = "Hide source code"

            title = nodes.paragraph(text=f"Function: {name}", classes=["toggle-button-title"])
            literal = nodes.literal_block(text=src_without_doc, language="python", classes=["toggle-content"])

            container += title
            container += literal
            output_nodes.append(container)

        return output_nodes

# ----------------------------
# InterlacedModuleDirective
# ----------------------------

class InterlacedModuleDirective(SphinxDirective):
    required_arguments = 1  # module path

    def run(self):
        module_path = self.arguments[0]
        module = importlib.import_module(module_path)
        output_nodes = []

        for name, func in inspect.getmembers(module, inspect.isfunction):
            # --- Render docstring via autodoc ---
            doc_nodes = FunctionDocumenter(self.env, self.name, self.options, self.state, self.lineno).get_doc(func)
            output_nodes.extend(doc_nodes)

            # --- Render source toggle ---
            src = inspect.getsource(func)
            container = nodes.container(classes=["toggle", "toggle-button"])
            container['data-toggle-label-show'] = "Show source code"
            container['data-toggle-label-hide'] = "Hide source code"
            literal = nodes.literal_block(text=src, language="python", classes=["toggle-content"])
            container += literal
            output_nodes.append(container)

        return output_nodes


def setup(app):
    app.add_directive("function_body", FunctionBodyDirective)
    app.add_directive("all_function_bodies", AllFunctionsDirective)
    app.add_directive("interlaced_module", InterlacedModuleDirective)
    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
