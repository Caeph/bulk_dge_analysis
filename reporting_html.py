from dominate import document
from dominate.tags import *
import io, base64, mimetypes
from dominate.util import raw
import os


class HTMLReporter:
    def __init__(self, title="Report"):
        self.doc = document(title=title)
        with self.doc.head:
            style(
                """
                body { margin: 0; padding: 20px; font-family: sans-serif; }
                .figure img { max-width: 100%; height: auto; display: block; margin: 10px 0; }
                pre { background-color: #f5f5f5; padding: 10px; border-radius: 4px; overflow-x: auto; }
                table { width: 100%; border-collapse: collapse; margin: 10px 0; }
                table th, table td { border: 1px solid #ddd; padding: 8px; text-align: left; }
                p.caption { font-size: 0.9em; color: #555; margin: 8px 0 4px; }
                """
            )
        with self.doc:
            h1(title)

    def start_section(self, section_title):
        with self.doc:
            h2(section_title)

    def send_figure(self, fig, caption, save_figure_path=None):
        if save_figure_path:
            os.makedirs(os.path.dirname(save_figure_path) or ".", exist_ok=True)
            fig.savefig(save_figure_path, bbox_inches="tight")

        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")

        data = base64.b64encode(buf.getvalue()).decode("ascii")
        with self.doc:
            with div():
                p(caption)
                img(src=f"data:image/png;base64,{data}")

    def send_text(self, text, caption=None, code_block=True):
        """
        Add plain text/log output.
        """
        with self.doc:
            if caption:
                p(caption, _class="caption")
            if code_block:
                pre(code(text))
            else:
                p(text)

    def send_image_from_path(self, image_path, caption):
        """
        Add a caption followed by an existing image referenced by file path.
        """
        with self.doc:
            # Caption first
            p(caption)
            # Then image by path
            with open(image_path, 'rb') as f:
                data = base64.b64encode(f.read()).decode('ascii')
            mime, _ = mimetypes.guess_type(image_path)
            src = f"data:{mime or 'image/png'};base64,{data}"
            div(
                img(src=src),
                _class="figure"
            )

    def save_table(self, table, coloring_rule, caption):
        with self.doc:
            p(caption)
            raw(table.style.map(coloring_rule).to_html(float_format='%.3E'))

    def save(self, filename="report.html"):
        with open(filename, "w") as f:
            f.write(self.doc.render())