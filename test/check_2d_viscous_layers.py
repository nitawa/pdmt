#!/usr/bin/env python3
"""Topology and viscous layer checks for PDMT's 2D VTU output."""

import argparse
from collections import Counter
import sys
import xml.etree.ElementTree as ET

VTK_LINE = 3
VTK_POLYGON = 7

def data_arrays(parent):
    return {
        array.get("Name"): (array.text or "").split()
        for array in parent.findall("./DataArray")
    }

def read_mesh(path):
    piece = ET.parse(path).getroot().find(".//Piece")
    if piece is None:
        raise ValueError("VTU has no UnstructuredGrid Piece")

    point_values = list(
        map(float, (piece.find("./Points/DataArray").text or "").split())
    )
    points = [
        tuple(point_values[index : index + 3])
        for index in range(0, len(point_values), 3)
    ]

    arrays = data_arrays(piece.find("./Cells"))
    connectivity = list(map(int, arrays["connectivity"]))
    offsets = list(map(int, arrays["offsets"]))
    types = list(map(int, arrays["types"]))

    cells = []
    begin = 0
    for end, cell_type in zip(offsets, types):
        cells.append((cell_type, connectivity[begin:end]))
        begin = end
    return points, cells

def check_viscous_layer(path, min_quad_count=1):
    points, cells = read_mesh(path)
    polygons = [nodes for cell_type, nodes in cells if cell_type == VTK_POLYGON]
    lines = [nodes for cell_type, nodes in cells if cell_type == VTK_LINE]

    if not polygons:
        raise AssertionError(f"{path}: no polygons found")

    quad_count = 0
    for cell_index, polygon in enumerate(polygons):
        if len(polygon) < 3 or len(set(polygon)) != len(polygon):
            raise AssertionError(f"{path}: polygon {cell_index} is degenerate")
        twice_area = 0.0
        for index, node in enumerate(polygon):
            next_node = polygon[(index + 1) % len(polygon)]
            x0, y0, _ = points[node]
            x1, y1, _ = points[next_node]
            twice_area += x0 * y1 - x1 * y0
        if abs(twice_area) < 1.0e-12:
            raise AssertionError(f"{path}: polygon {cell_index} has zero area")
        if len(polygon) == 4:
            quad_count += 1

    if quad_count < min_quad_count:
        raise AssertionError(
            f"{path}: expected at least {min_quad_count} quad cells for viscous layers, found {quad_count}"
        )
    print(f"PASS: {path} verified with {len(polygons)} polygons ({quad_count} quad cells).")

def main():
    parser = argparse.ArgumentParser(description="Check 2D dual mesh with viscous layers.")
    parser.add_argument("vtu_file", help="Path to output VTU file")
    parser.add_argument("--min-quads", type=int, default=1, help="Minimum expected quad cells")
    args = parser.parse_args()
    check_viscous_layer(args.vtu_file, min_quad_count=args.min_quads)

if __name__ == "__main__":
    main()
