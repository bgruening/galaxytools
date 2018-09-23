import sys
import argparse
import plotly
import plotly.graph_objs as go
import pandas as pd

def main(infile, col_dimensions, col_color):
    """
    Produce an interactive paracords plotting html
    Args:
        infile: str, tabular file
        col_dimensions: str, comma separated index numbers. For example: "3,4,5"
        col_color: str, index number
    """
    df = pd.read_csv(infile, sep='\t', parse_dates=True)

    dimensions = []
    col_dimensions = [int(x)-1 for x in col_dimensions.split(',')]
    for col in col_dimensions:
        values = df[df.columns[col]]
        if all(type(e) is int for e in values ):
            dimensions.append(
                dict(   values = values,
                        tickformat = ",.2r",
                        label = df.columns[col])
            )
        elif all(type(e) is float for e in values ):
            dimensions.append(
                dict(   values = values,
                        tickformat = "g",
                        label = df.columns[col])
            )
        else:
            unique_values = list(set(values))
            dimensions.append(
                dict(   range = [0, len(unique_values)-1],
                        tickvals = list(range(len(unique_values))),
                        ticktext = [str(e) for e in unique_values],
                        values = list(map(lambda e: unique_values.index(e), values )),
                        label = df.columns[col])
            )

    col_color = int(col_color) - 1
    colors = df[df.columns[col_color]]
    if all(type(e) is int for e in colors ):
        tickformat = ",.2r"
    elif all(type(e) is float for e in colors ):
        tickformat = "g"
    else:
        sys.exit("Error: the column for coloring must contain all numerical values!")

    dimensions.append(
        dict(
                values = colors,
                tickformat = tickformat,
                label = df.columns[col_color]
        )
    )

    line = dict(
                color = colors,
                colorscale = 'Jet',
                showscale = True,
                reversescale = True
    )

    data = [
            go.Parcoords(
                line = line,
                dimensions = dimensions
            )
    ]

    plotly.offline.plot(data, filename = "output.html", auto_open=False)

if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument( "-i", "--input", dest="infile", required=True)
    aparser.add_argument( "-d", "--col_dimensions", dest="col_dimensions")
    aparser.add_argument( "-c", "--col_color", dest="col_color")
    args = aparser.parse_args()

    main(args.infile, args.col_dimensions, args.col_color)