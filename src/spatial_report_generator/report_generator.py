"""
Python class containing methods to populate report template.
"""
import base64
from io import BytesIO

import jinja2
import matplotlib.pyplot as plt
import plotly
import seaborn as sns

from src.utils.logger import Logger

logger = Logger()


class ReportGenerator:
    """Python class to populate report template with generated figures."""

    def __init__(self, template_path):
        """Initializing ReportGenerator instance variables.

        Args:
            template_path (str): path to report template in cache.
        """
        self.template_path = template_path
        self.environment = jinja2.Environment()  # NOSONAR
        with open(template_path, "r", encoding="utf-8") as f:
            self.report = self.environment.from_string(f.read())

        f.close()
        logger.info(f"Report generator initialized")

    def fig_to_base64(self, fig):
        """Convert matplotlib figure to base64 string

        Args:
            fig (matplotlib.figure.Figure): matplotlib figure

        Returns:
            str: base64 string
        """

        buffer = BytesIO()
        fig.savefig(buffer, format="png")
        buffer.seek(0)
        return base64.b64encode(buffer.read()).decode("utf-8")

    def format_to_html(self, fig):
        if fig is not None:
            # Write the figure
            if isinstance(fig, plt.Figure):
                html_img = (
                    f'<img src="data:image/png;base64,{self.fig_to_base64(fig)}">'
                )
            elif isinstance(fig, plotly.graph_objs._figure.Figure):
                html_img = plotly.io.to_html(
                    fig, include_plotlyjs=False, include_mathjax="cdn", full_html=False
                )
            else:
                html_img = fig
            return html_img
        return ""

    def populate_report(self, **kwargs):
        """Populates the report template using the key word arguments defined. Takes as input the
        parameters to populate the template with.

        Returns:
            str: populated report html.
        """
        logger.info("<populate_report> Populating report")
        return self.report.render(kwargs)
