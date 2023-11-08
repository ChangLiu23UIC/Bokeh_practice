from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, show, output_file

"""
Plot a protein coverage map with the input uniport id 
"""

def protein_coverage_map(uniport: str, found_dict: dict, protein_list: list, uniport_list: list):
    # find the protein sequence and the psm first with uniport ID
    protein_seq = protein_list[uniport_list.index(uniport)]
    psm_list = found_dict[uniport]

    # generate a plot in html
    output_file("Result.html")

    # get the position of psm
    peptide_positions = [(protein_seq.index(pep), protein_seq.index(pep) + len(pep)) for pep in psm_list]

    # Create a ColumnDataSource with the peptide positions and names
    source = ColumnDataSource(data = dict(
        # start pos
        left = [pos[0] for pos in peptide_positions],
        # end pos
        right = [pos[1] for pos in peptide_positions],
        # name of the psm
        name = psm_list,
        # position of psm in graph on y axis
        y_pos = [i + 1 for i in range(len(psm_list))]
    ))

    # Create a figure as a background to put plot on it
    plot = figure(
        x_range=(0, len(protein_seq)),
        y_range=(0, len(psm_list)),
        title = "Protein Coverage Map",
        x_axis_label = 'Protein sequence position',
        y_axis_label = 'PSM',
        width = 1080,
        height = 1000)

    # Add interval bars
    plot.hbar(y = 'y_pos', left = 'left', right = 'right', height=0.8, source=source, fill_color='green', line_color='white')

    # change the label index (one at a time) to show every psm in y axis
    plot.yaxis.ticker = [i+1 for i in range(len(psm_list))]
    plot.yaxis.major_label_overrides = {i+1: name for i, name in enumerate(psm_list)}
    plot.yaxis.ticker = [i+1 for i in range(len(protein_seq))]
    plot.xaxis.major_label_overrides = {i+1: name for i, name in enumerate(protein_seq)}

    show(plot)



