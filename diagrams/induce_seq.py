import os
import textwrap

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUTPUT_JPEG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "induce_seq_diagrams.jpg")

fig, ax = plt.subplots()

STRAND_WIDTH = 0.06
STRAND_SEPARATION = 0.04
STRAND_LENGTH = 1.3
STRAND_OVERHANG = 0.4
DSB_WIDTH = 0.05
LABEL_OFFSET = 0.05
LABEL_VERTICAL_OFFSET = 0.05
BREAK_POSITIONS = [ 0.1, 0.95]
STRAND_EDGE_WIDTH = 1
LABEL_FONT_SIZE = 15
# Font size for labels that aren't attached to an arrow (e.g. strand-end 5'/3' labels, the DSB
# label) -- kept separate from LABEL_FONT_SIZE so the two can be sized independently.
NON_ARROW_LABEL_FONT_SIZE = 12
DENATURED_STRANDS_HORIZONTAL_OFFSET = -0.7
DENATURED_STRANDS_VERTICAL_OFFSET = 0.3

OTHER_FRAGMENT_LENGTHS = [0.3, 0.6, 0.9]
FRAGMENT_ROW_SPACING = 0.25

# Vertical shift applied to everything below diagram 4 -- diagrams 5, 6, 7 (the fragments pair
# and the flow cell) and all of their arrows/labels -- as a single block. fragments_position_y
# (below) anchors that whole block: flow_cell_position_y and every arrow/label positioned below
# diagram 4 are computed relative to it, so adding this here shifts the whole block up (positive)
# or down (negative) as one unit without disturbing the gaps within it.
BELOW_DIAGRAM_4_VERTICAL_SHIFT = 0.2

# Labels for the arrows drawn between consecutive diagrams in __main__, top-to-bottom. Edit
# freely -- there's one arrow per gap between diagrams, so this needs exactly 5 entries (the
# draw_dna_fragments / draw_dna_fragments_with_p7_adapters pair counts as one, since they're
# drawn side by side rather than as separate rows -- see PAIR_ARROW_LABEL for the sideways arrow
# between those two specifically).
ARROW_LABELS = [
    "DSB end blunting",
    "P5 adapter ligation",
    "fragmentation",
    "fragmentation cont.",
    "denaturation and binding of fragments to a flow cell",
]

# Label for the sideways arrow between draw_dna_fragments and draw_dna_fragments_with_p7_adapters.
PAIR_ARROW_LABEL = "P7 adapter ligation and size filtering"
# Label pointing out the tail (dark magenta + dark orange) P7 adapter on diagram 6's first
# fragment, placed below PAIR_ARROW_LABEL.
HALF_FUNCTIONAL_P7_LABEL = "half-functional P7 adapter"
# How far above the diagrams 5/6 pair's shared center y the arrow between them (and its label)
# sits -- kept separate from fragments_position_y so raising the arrow doesn't also raise the
# diagrams themselves.
PAIR_ARROW_VERTICAL_OFFSET = 0.2


def draw_rectangle(x, y, width, height, **kwargs):
    """Draws a rectangle with its bottom-left corner at (x, y) on the module-level figure."""
    rect = Rectangle((x, y), width, height, **kwargs)
    ax.add_patch(rect)
    return rect


def draw_strand_segment(x, y, width, height, skip_edges=(), facecolor="lightgray", edgecolor="black", linewidth=STRAND_EDGE_WIDTH):
    """
    Draws a rectangle representing a strand segment: filled with facecolor, with a border of
    edgecolor on every edge except those named in skip_edges (any of 'left', 'right', 'top',
    'bottom'). Rectangle's own edgecolor applies uniformly to all four sides, so to leave some
    edges open the fill is drawn borderless and each wanted edge is drawn on top as its own line.
    """
    rect = draw_rectangle(x, y, width, height, facecolor=facecolor, edgecolor="none")
    edges = {
        "bottom": ((x, x + width), (y, y)),
        "top": ((x, x + width), (y + height, y + height)),
        "left": ((x, x), (y, y + height)),
        "right": ((x + width, x + width), (y, y + height)),
    }
    for edge, (xs, ys) in edges.items():
        if edge not in skip_edges:
            ax.plot(xs, ys, color=edgecolor, linewidth=linewidth, solid_capstyle="projecting")
    return rect


def draw_end_label(x, y, label, side):
    """
    Draws a strand-end label (e.g. "5'" or "3'") vertically centered at y, offset LABEL_OFFSET
    to the given side ('left' or 'right') of x, so it sits just outside the strand end at x.
    """
    ha = "right" if side == "left" else "left"
    offset = -LABEL_OFFSET if side == "left" else LABEL_OFFSET
    ax.text(x + offset, y, label, ha=ha, va="center", fontsize=NON_ARROW_LABEL_FONT_SIZE)


ARROW_LENGTH = 0.2
LAST_ARROW_LENGTH = 0.35
ARROW_LABEL_OFFSET = 0.1


def draw_labeled_arrow(x, y_center, label, length=ARROW_LENGTH, label_side="right", label_dx=0, label_dy=0):
    """
    Draws a short downward-pointing arrow centered at (x, y_center), with label as text just to
    its side (label_side, "left" or "right"), vertically centered on the arrow. label_dx/label_dy
    additionally shift just the label (not the arrow) from that default position, for callers
    that need the caption's position independently adjustable from the arrow's. Used to mark the
    transition between consecutive diagrams in __main__.
    """
    ax.annotate(
        "", xy=(x, y_center - length / 2), xytext=(x, y_center + length / 2),
        arrowprops=dict(arrowstyle="->", color="black", linewidth=1.5),
    )
    ha = "left" if label_side == "right" else "right"
    offset = ARROW_LABEL_OFFSET if label_side == "right" else -ARROW_LABEL_OFFSET
    ax.text(x + offset + label_dx, y_center + label_dy, label, ha=ha, va="center", fontsize=LABEL_FONT_SIZE)


def draw_horizontal_labeled_arrow(x_center, y, label=None, length=ARROW_LENGTH):
    """
    Draws a short rightward-pointing arrow centered at (x_center, y), with label as text just
    above it, horizontally centered on the arrow, if given. The horizontal counterpart to
    draw_labeled_arrow, used to mark the transition between two side-by-side diagrams.
    """
    ax.annotate(
        "", xy=(x_center + length / 2, y), xytext=(x_center - length / 2, y),
        arrowprops=dict(arrowstyle="->", color="black", linewidth=1.5),
    )
    if label:
        ax.text(x_center, y + ARROW_LABEL_OFFSET, label, ha="center", va="bottom", fontsize=LABEL_FONT_SIZE)


def draw_two_tone_strand(x, y, length, height, first_color, second_color):
    """
    Draws a horizontal strand-shaped rectangle of the given total length as two equal-length
    halves side by side -- first_color on the left half, second_color on the right half --
    with no border on the seam between them (the same physical piece), but a normal border
    everywhere else (top, bottom, and both outer ends).
    """
    half_length = length / 2
    draw_strand_segment(x, y, half_length, height, skip_edges=("right",), facecolor=first_color)
    draw_strand_segment(x + half_length, y, half_length, height, skip_edges=("left",), facecolor=second_color)


def draw_dsb(position=(0, 0)):
    """
    Draws a double strand break centered at position=(x, y): two strand segments (top and
    bottom, STRAND_WIDTH thick, separated by STRAND_SEPARATION) leading up to the break from
    the left, and two more continuing on from the break to the right. Every coordinate is
    computed relative to position, so the whole diagram can be moved around just by changing it.

    The left two segments share a left edge (left open/no border, since the strand continues
    off frame there), but the top one's break end sticks out a bit further right than the
    bottom one's (by STRAND_OVERHANG), showing a staggered break. The right two segments each
    start DSB_WIDTH to the right of their matching left segment's break end, and share a right
    edge (also open/no border, continuing off frame) -- so the top-right segment ends up
    shorter than the bottom-right one, by that same stagger.
    """
    pos_x, pos_y = position

    bottom_y = pos_y - STRAND_SEPARATION / 2 - STRAND_WIDTH
    top_y = pos_y + STRAND_SEPARATION / 2

    left_x = pos_x - STRAND_LENGTH
    right_x = pos_x + STRAND_LENGTH

    bottom_break_x = pos_x
    top_break_x = bottom_break_x + STRAND_OVERHANG

    # Left segments: left edge open (shared, off frame); right edge is the break end.
    draw_strand_segment(left_x, bottom_y, bottom_break_x - left_x, STRAND_WIDTH, skip_edges=("left",))
    draw_strand_segment(left_x, top_y, top_break_x - left_x, STRAND_WIDTH, skip_edges=("left",))

    # Right segments: left edge is DSB_WIDTH to the right of the matching left segment's break
    # end; right edge open (shared, off frame).
    draw_strand_segment(bottom_break_x + DSB_WIDTH, bottom_y, right_x - (bottom_break_x + DSB_WIDTH), STRAND_WIDTH, skip_edges=("right",))
    draw_strand_segment(top_break_x + DSB_WIDTH, top_y, right_x - (top_break_x + DSB_WIDTH), STRAND_WIDTH, skip_edges=("right",))

    # Antiparallel strands: the top strand runs 5'->3' left to right, the bottom strand 3'->5'.
    draw_end_label(left_x, top_y + STRAND_WIDTH / 2 + LABEL_VERTICAL_OFFSET, "5'", side="left")
    draw_end_label(right_x, top_y + STRAND_WIDTH / 2 + LABEL_VERTICAL_OFFSET, "3'", side="right")
    draw_end_label(left_x, bottom_y + STRAND_WIDTH / 2 - LABEL_VERTICAL_OFFSET, "3'", side="left")
    draw_end_label(right_x, bottom_y + STRAND_WIDTH / 2 - LABEL_VERTICAL_OFFSET, "5'", side="right")

    ax.text(
        bottom_break_x + STRAND_OVERHANG / 2, top_y + STRAND_WIDTH / 2 + ARROW_LABEL_OFFSET,
        "double strand break (DSB)", ha="center", va="bottom", fontsize=NON_ARROW_LABEL_FONT_SIZE,
    )


def draw_blunted_dsb(position=(0, 0), fragment_shift=0):
    """
    Draws a blunted (no-overhang) double strand break centered at position=(x, y): like
    draw_dsb, but both breaks are flush instead of staggered, so each pair of parallel
    rectangles ends up with aligned left and right edges. Only one rectangle per pair is
    shortened to achieve this (moving just the one edge that needs to move) -- the other keeps
    the exact position and size it has in draw_dsb.

    fragment_shift, if nonzero, extends each piece's break end further out by that amount (the
    left piece's break end moves left, the right piece's moves right), widening the gap between
    the two DNA fragments. The far/outer ends (which continue off frame) stay put.
    """
    pos_x, pos_y = position

    bottom_y = pos_y - STRAND_SEPARATION / 2 - STRAND_WIDTH
    top_y = pos_y + STRAND_SEPARATION / 2

    left_x = pos_x - STRAND_LENGTH
    right_x = pos_x + STRAND_LENGTH

    break_x = pos_x - fragment_shift
    right_segment_start_x = pos_x + STRAND_OVERHANG + DSB_WIDTH + fragment_shift

    # Left segments: left edge open (shared, off frame); right edge is the break end. Bottom-left
    # is unmodified from draw_dsb; top-left is shortened (its right edge moved left) to match it.
    draw_strand_segment(left_x, bottom_y, break_x - left_x, STRAND_WIDTH, skip_edges=("left",))
    draw_strand_segment(left_x, top_y, break_x - left_x, STRAND_WIDTH, skip_edges=("left",))

    # Right segments: left edge is the break end. Top-right is unmodified from draw_dsb;
    # bottom-right is shortened (its left edge moved right) to match it. Right edge open
    # (shared, off frame).
    draw_strand_segment(right_segment_start_x, bottom_y, right_x - right_segment_start_x, STRAND_WIDTH, skip_edges=("right",))
    draw_strand_segment(right_segment_start_x, top_y, right_x - right_segment_start_x, STRAND_WIDTH, skip_edges=("right",))

    # Antiparallel strands: the top strand runs 5'->3' left to right, the bottom strand 3'->5'.
    # (Not labeled here -- only draw_dsb, the first diagram, gets 5'/3' text labels.)


ADAPTER_COLOR_ORANGE_DARK = "#B35900"
ADAPTER_COLOR_ORANGE_LIGHT = "#FFCC80"
ADAPTER_COLOR_MAGENTA_DARK = "#8B008B"
ADAPTER_COLOR_MAGENTA_LIGHT = "#FFB3E6"
ADAPTER_COLOR_BLUE_DARK = "#003366"
ADAPTER_COLOR_BLUE_LIGHT = "#99CCFF"
ADAPTER_COLOR_GREEN_DARK = "#1B5E20"
ADAPTER_COLOR_GREEN_LIGHT = "#B2FFB2"
ADAPTER_SET_GAP = 0.1
FRAGMENT_SHIFT = 0.1


def draw_blunted_dsb_with_adapters(position=(0, 0)):
    """
    Draws a draw_blunted_dsb diagram (with its two DNA fragments shifted apart by
    FRAGMENT_SHIFT, to leave room for longer adapters) with a P5 adapter ligated to each of its
    two blunted DNA ends (the left piece's right end, and the right piece's left end, i.e. the
    two sides of the break). Each ligated adapter is itself double-stranded, so it's drawn as
    two rectangles the same shape as a DNA strand, butted up against the strand end they're
    ligated to -- unlike the open/borderless strand ends used elsewhere in these diagrams for
    "continues off frame", the adapter-DNA junction keeps its border, since this is a real,
    drawn boundary between two distinct, fully-drawn pieces. Each of the two adapter strands
    (top and bottom) is itself two equal-length sections, blue then green, with no border on
    the seam between them (the same physical piece); the two adapters at each break are
    centered in what was the break gap, leaving ADAPTER_SET_GAP of open space between them.
    Whichever strand's cut end is 5' at that junction uses the dark blue/green shades, and
    whichever is 3' uses the light shades, so the two opposite strands read as obviously
    different at a glance.
    """
    draw_blunted_dsb(position, fragment_shift=FRAGMENT_SHIFT)

    pos_x, pos_y = position
    bottom_y = pos_y - STRAND_SEPARATION / 2 - STRAND_WIDTH
    top_y = pos_y + STRAND_SEPARATION / 2

    break_x = pos_x - FRAGMENT_SHIFT
    right_segment_start_x = pos_x + STRAND_OVERHANG + DSB_WIDTH + FRAGMENT_SHIFT
    adapter_length = (right_segment_start_x - break_x - ADAPTER_SET_GAP) / 2

    # Left DNA end: the top strand's cut end is 3' (light), the bottom strand's is 5' (dark).
    draw_two_tone_strand(break_x, top_y, adapter_length, STRAND_WIDTH, ADAPTER_COLOR_BLUE_LIGHT, ADAPTER_COLOR_GREEN_LIGHT)
    draw_two_tone_strand(break_x, bottom_y, adapter_length, STRAND_WIDTH, ADAPTER_COLOR_BLUE_DARK, ADAPTER_COLOR_GREEN_DARK)

    # Right DNA end: the top strand's cut end is 5' (dark), the bottom strand's is 3' (light).
    draw_two_tone_strand(right_segment_start_x - adapter_length, top_y, adapter_length, STRAND_WIDTH, ADAPTER_COLOR_BLUE_DARK, ADAPTER_COLOR_GREEN_DARK)
    draw_two_tone_strand(right_segment_start_x - adapter_length, bottom_y, adapter_length, STRAND_WIDTH, ADAPTER_COLOR_BLUE_LIGHT, ADAPTER_COLOR_GREEN_LIGHT)


def split_spans_at_breaks(spans, break_xs, gap):
    """
    spans is an ordered, non-overlapping list of (x_start, x_end, facecolor, seamless_before)
    tuples (which may already have gaps, or borderless seams marked by seamless_before, between
    some of them -- see draw_strand_row). For each x in break_xs that falls strictly inside one
    of the spans, splits that span into two, leaving `gap` of empty space centered on x (each
    side shortened by gap / 2); the second half is never seamless with what preceded it, since
    there's now a real gap there. Returns the resulting list of spans, in the same order.
    """
    for break_x in break_xs:
        new_spans = []
        for x_start, x_end, facecolor, seamless_before in spans:
            if x_start < break_x < x_end:
                new_spans.append((x_start, break_x - gap / 2, facecolor, seamless_before))
                new_spans.append((break_x + gap / 2, x_end, facecolor, False))
            else:
                new_spans.append((x_start, x_end, facecolor, seamless_before))
        spans = new_spans
    return spans


def draw_strand_row(y, spans):
    """
    Draws a horizontal strand row from an ordered list of (x_start, x_end, facecolor,
    seamless_before) spans: every span is fully bordered (a real break/junction) except the
    left edge of the first span and the right edge of the last span, which stay open,
    continuing off frame (matching the open outer ends used everywhere else in these diagrams),
    and except a seam where a span's seamless_before is True, which leaves both sides of that
    seam borderless (the same physical piece, e.g. the two halves of a two-tone adapter section).
    """
    for i, (x_start, x_end, facecolor, seamless_before) in enumerate(spans):
        skip_edges = []
        if i == 0 or seamless_before:
            skip_edges.append("left")
        if i == len(spans) - 1 or (i + 1 < len(spans) and spans[i + 1][3]):
            skip_edges.append("right")
        draw_strand_segment(x_start, y, x_end - x_start, STRAND_WIDTH, skip_edges=tuple(skip_edges), facecolor=facecolor)


def draw_fragment_dna(position=(0, 0), break_positions=()):
    """
    Draws the same picture as draw_blunted_dsb_with_adapters (a blunted DSB, its two DNA
    fragments shifted apart by FRAGMENT_SHIFT, with a blue/green P5 adapter ligated at both
    ends), but with the DNA additionally fragmented at break_positions: each entry is a fraction from 0
    (the diagram's left end) to 1 (its right end) of the whole strand's length, giving the x
    position of a break. Every break cuts straight across both the top and bottom strand at the
    same x (unlike the original DSB break, which is staggered) -- random-shearing fragmentation
    breaks both strands cleanly at the same point.
    """
    pos_x, pos_y = position
    bottom_y = pos_y - STRAND_SEPARATION / 2 - STRAND_WIDTH
    top_y = pos_y + STRAND_SEPARATION / 2

    left_x = pos_x - STRAND_LENGTH
    right_x = pos_x + STRAND_LENGTH

    break_x = pos_x - FRAGMENT_SHIFT
    right_segment_start_x = pos_x + STRAND_OVERHANG + DSB_WIDTH + FRAGMENT_SHIFT
    adapter_length = (right_segment_start_x - break_x - ADAPTER_SET_GAP) / 2
    half_adapter_length = adapter_length / 2

    # Each P5 adapter is drawn as two equal-length, seamlessly-joined halves: blue then green.
    # The top strand's cut end is 3' (light) on the left, 5' (dark) on the right; the bottom
    # strand is the opposite -- 5' (dark) on the left, 3' (light) on the right.
    top_spans = [
        (left_x, break_x, "lightgray", False),
        (break_x, break_x + half_adapter_length, ADAPTER_COLOR_BLUE_LIGHT, False),
        (break_x + half_adapter_length, break_x + adapter_length, ADAPTER_COLOR_GREEN_LIGHT, True),
        (right_segment_start_x - adapter_length, right_segment_start_x - half_adapter_length, ADAPTER_COLOR_BLUE_DARK, False),
        (right_segment_start_x - half_adapter_length, right_segment_start_x, ADAPTER_COLOR_GREEN_DARK, True),
        (right_segment_start_x, right_x, "lightgray", False),
    ]
    bottom_spans = [
        (left_x, break_x, "lightgray", False),
        (break_x, break_x + half_adapter_length, ADAPTER_COLOR_BLUE_DARK, False),
        (break_x + half_adapter_length, break_x + adapter_length, ADAPTER_COLOR_GREEN_DARK, True),
        (right_segment_start_x - adapter_length, right_segment_start_x - half_adapter_length, ADAPTER_COLOR_BLUE_LIGHT, False),
        (right_segment_start_x - half_adapter_length, right_segment_start_x, ADAPTER_COLOR_GREEN_LIGHT, True),
        (right_segment_start_x, right_x, "lightgray", False),
    ]

    break_xs = [left_x + fraction * (right_x - left_x) for fraction in break_positions]

    draw_strand_row(top_y, split_spans_at_breaks(top_spans, break_xs, DSB_WIDTH))
    draw_strand_row(bottom_y, split_spans_at_breaks(bottom_spans, break_xs, DSB_WIDTH))

    # Antiparallel strands: the top strand runs 5'->3' left to right, the bottom strand 3'->5'.
    # (Not labeled here -- only draw_dsb, the first diagram, gets 5'/3' text labels.)



def draw_dna_fragment_row(y, x_start, spans):
    """
    Draws one row (one strand) of a standalone, closed DNA fragment at height y, starting at
    x_start, from an ordered list of (length, facecolor, seamless_before) segments. Every
    segment is fully bordered, except where a segment's seamless_before is True, which leaves
    both sides of that seam borderless (the same physical piece, e.g. the two halves of a
    two-tone adapter section) -- unlike draw_strand_row, no segment's outer edge is left open
    for "continues off frame", since this represents an already-isolated, complete fragment.
    """
    x = x_start
    for i, (length, facecolor, seamless_before) in enumerate(spans):
        skip_edges = []
        if seamless_before:
            skip_edges.append("left")
        if i + 1 < len(spans) and spans[i + 1][2]:
            skip_edges.append("right")
        draw_strand_segment(x, y, length, STRAND_WIDTH, skip_edges=tuple(skip_edges), facecolor=facecolor)
        x += length


def fragmented_dna_lengths():
    """
    Returns (left_dna_length, right_dna_length): the lengths of the DNA piece immediately
    adjacent to each P5 adapter after fragmentation at BREAK_POSITIONS -- i.e. the same size
    the equivalent gray segment ends up in draw_fragment_dna (from the nearest break to the
    adapter's attachment point), rather than the full, unfragmented piece draw_blunted_dsb and
    draw_blunted_dsb_with_adapters use.
    """
    left_x = -STRAND_LENGTH
    right_x = STRAND_LENGTH
    break_x = -FRAGMENT_SHIFT
    right_segment_start_x = STRAND_OVERHANG + DSB_WIDTH + FRAGMENT_SHIFT

    break_xs = sorted(left_x + fraction * (right_x - left_x) for fraction in BREAK_POSITIONS)

    left_breaks = [x for x in break_xs if x < break_x]
    left_boundary = left_breaks[-1] + DSB_WIDTH / 2 if left_breaks else left_x
    left_dna_length = break_x - left_boundary

    right_breaks = [x for x in break_xs if x > right_segment_start_x]
    right_boundary = right_breaks[0] - DSB_WIDTH / 2 if right_breaks else right_x
    right_dna_length = right_boundary - right_segment_start_x

    return left_dna_length, right_dna_length


def build_dna_fragments():
    """
    Builds the list of fragments used by draw_dna_fragments: each entry is (top spans, bottom
    spans, label its left end, label its right end), with spans as ordered (length, facecolor,
    seamless_before) triples, left to right (see draw_dna_fragment_row). The two adapter-bearing
    fragments only get a 5'/3' text label on their plain-DNA end -- their adapter end is already
    identified by its color, matching how draw_fragment_dna doesn't label the adapters either.
    Each P5 adapter is two equal-length, seamlessly-joined halves: blue then green. Their DNA
    lengths match the fragmented (post-break) pieces from draw_fragment_dna (see
    fragmented_dna_lengths), not the full, unfragmented piece. Factored out so other diagrams
    (e.g. draw_dna_fragments_with_p7_adapters) can reuse the exact same fragments, lengths, and
    orientations.
    """
    adapter_length = (STRAND_OVERHANG + DSB_WIDTH + 2 * FRAGMENT_SHIFT - ADAPTER_SET_GAP) / 2
    half_adapter_length = adapter_length / 2
    left_dna_length, right_dna_length = fragmented_dna_lengths()

    return [
        (
            [(left_dna_length, "lightgray", False), (half_adapter_length, ADAPTER_COLOR_BLUE_LIGHT, False), (half_adapter_length, ADAPTER_COLOR_GREEN_LIGHT, True)],
            [(left_dna_length, "lightgray", False), (half_adapter_length, ADAPTER_COLOR_BLUE_DARK, False), (half_adapter_length, ADAPTER_COLOR_GREEN_DARK, True)],
            True, False,
        ),
        (
            [(half_adapter_length, ADAPTER_COLOR_BLUE_DARK, False), (half_adapter_length, ADAPTER_COLOR_GREEN_DARK, True), (right_dna_length, "lightgray", False)],
            [(half_adapter_length, ADAPTER_COLOR_BLUE_LIGHT, False), (half_adapter_length, ADAPTER_COLOR_GREEN_LIGHT, True), (right_dna_length, "lightgray", False)],
            False, True,
        ),
    ] + [
        ([(length, "lightgray", False)], [(length, "lightgray", False)], True, True)
        for length in OTHER_FRAGMENT_LENGTHS
    ]


def dna_fragments_max_width():
    """Returns the length of the longest fragment build_dna_fragments produces."""
    return max(sum(length for length, _, _ in top_spans) for top_spans, _, _, _ in build_dna_fragments())


def dna_fragments_row_width():
    """
    Returns the total width of the row of fragment columns draw_dna_fragments lays out -- i.e.
    from the first fragment's column to the last one's.
    """
    return (len(build_dna_fragments()) - 1) * FRAGMENT_ROW_SPACING


def dna_fragments_with_p7_row_width():
    """
    Returns the total width of the row of fragment columns draw_dna_fragments_with_p7_adapters
    lays out. Narrower than dna_fragments_row_width() by one FRAGMENT_ROW_SPACING per excluded
    fragment, since (unlike that function) the remaining columns close up the gap a skipped
    fragment would otherwise leave, rather than keeping their original, wider-spaced positions.
    """
    return (len(build_dna_fragments()) - len(EXCLUDED_FRAGMENT_INDICES) - 1) * FRAGMENT_ROW_SPACING


def draw_dna_fragments(position=(0, 0)):
    """
    Draws a row of standalone, complete DNA fragments, each column to the right of the previous
    one (spaced FRAGMENT_ROW_SPACING apart), the whole row horizontally centered on position's
    x and each column vertically centered on position's y. Unlike the other diagrams in this
    file, every fragment here is fully bordered on all sides: each one is a separate, complete
    molecule, rather than a piece of a longer sequence continuing off frame.

    Two of the fragments have the same lengths and orientations as the two DNA pieces (each
    with one ligated P5 adapter) from draw_fragment_dna / draw_blunted_dsb_with_adapters -- just
    not at the same positions, since here they're independent columns in this row. The
    remaining fragments are plain (no adapter), with lengths taken from the
    OTHER_FRAGMENT_LENGTHS global.
    """
    pos_x, pos_y = position
    row_left_x = pos_x - dna_fragments_row_width() / 2

    for i, (top_spans, bottom_spans, label_left, label_right) in enumerate(build_dna_fragments()):
        column_x = row_left_x + i * FRAGMENT_ROW_SPACING
        right_x = column_x + STRAND_SEPARATION / 2 + STRAND_WIDTH / 2
        left_x = column_x - STRAND_SEPARATION / 2 - STRAND_WIDTH / 2
        fragment_length = sum(length for length, _, _ in top_spans)
        bottom_y = pos_y - fragment_length / 2

        # Antiparallel strands: the top strand runs 5'->3' bottom to top, the bottom 3'->5' --
        # drawn as the left and right columns respectively (mirrored from their build_dna_fragments
        # order so each fragment's light/near-3' adapter half ends up on the left; not labeled
        # here -- only draw_dsb, the first diagram, gets 5'/3' text labels).
        draw_dna_fragment_column(left_x, bottom_y, top_spans)
        draw_dna_fragment_column(right_x, bottom_y, bottom_spans)


P7_ADAPTER_LENGTH = 0.3
EXCLUDED_FRAGMENT_INDICES = [2]


def p7_adapter_extension(is_tail, reversed_order=False):
    """
    Returns the (length, facecolor, seamless_before) segments of a P7 adapter attached to one
    end of a strand, in near-DNA-to-far order: just a light magenta half if not is_tail, or a
    dark magenta half followed (seamlessly) by a dark orange half if is_tail. If reversed_order,
    returns them far-to-near instead (with the seamless_before flags recomputed to match), for
    prepending/extending before a DNA span rather than appending/extending after one.
    """
    half = P7_ADAPTER_LENGTH / 2
    if not is_tail:
        return [(half, ADAPTER_COLOR_MAGENTA_LIGHT, False)]
    near = (half, ADAPTER_COLOR_MAGENTA_DARK, False)
    far = (half, ADAPTER_COLOR_ORANGE_DARK, True)
    if reversed_order:
        return [(far[0], far[1], False), (near[0], near[1], True)]
    return [near, far]


def add_p7_adapters(spans, is_top, label_left, label_right):
    """
    Returns spans with a P7 adapter spliced onto whichever end(s) are plain DNA (label_left
    and/or label_right), using the antiparallel-strand rule for which end is the 5' (tail) one:
    at a left end the top strand is 5', at a right end the bottom strand is.
    """
    result = list(spans)
    if label_left:
        result = p7_adapter_extension(is_tail=is_top, reversed_order=True) + result
    if label_right:
        result = result + p7_adapter_extension(is_tail=not is_top, reversed_order=False)
    return result


def dna_fragments_with_p7_max_width():
    """Returns the length of the longest fragment draw_dna_fragments_with_p7_adapters produces."""
    return max(
        sum(length for length, _, _ in add_p7_adapters(top_spans, True, label_left, label_right))
        for top_spans, _, label_left, label_right in build_dna_fragments()
    )


def dna_fragments_pair_extents():
    """
    Returns (top_offset, bottom_offset): how far above and below their shared center y the
    columns drawn by draw_dna_fragments and draw_dna_fragments_with_p7_adapters actually reach.

    Unlike dna_fragments_max_width() / dna_fragments_with_p7_max_width(), which give each
    fragment's total column length, this isn't just half of that: draw_dna_fragments_with_p7_adapters
    centers every column on the DNA-only fragment_length (see its bottom_y/top_y), then extends a
    P7 adapter out from just one end for the two adapter-bearing fragments -- so those columns
    reach further to one side of center than the other, and the largest such one-sided reach can
    exceed half of the largest fragment's total length. Callers that need to bound or place
    something (an arrow, a label, the next diagram down) relative to this pair should use these
    exact offsets rather than assuming symmetry.
    """
    top_offset = dna_fragments_max_width() / 2
    bottom_offset = dna_fragments_max_width() / 2
    for i, (top_spans, _, label_left, label_right) in enumerate(build_dna_fragments()):
        if i in EXCLUDED_FRAGMENT_INDICES:
            continue
        fragment_length = sum(length for length, _, _ in top_spans)
        top_offset = max(top_offset, fragment_length / 2 + (P7_ADAPTER_LENGTH if label_right else 0))
        bottom_offset = max(bottom_offset, fragment_length / 2 + (P7_ADAPTER_LENGTH if label_left else 0))
    return top_offset, bottom_offset


def draw_dna_fragments_with_p7_adapters(position=(0, 0)):
    """
    Draws the same fragment row as draw_dna_fragments, except:
    - any fragment index in EXCLUDED_FRAGMENT_INDICES is skipped entirely, and the remaining
      fragments close up the gap that would otherwise leave -- drawn one to a compacted column
      index rather than their original index i, over the (narrower) dna_fragments_with_p7_row_width().
    - every DNA end that doesn't already have a P5 adapter (i.e. every end draw_dna_fragments
      would otherwise leave as a plain, unlabeled end) is instead ligated with a P7 adapter (see
      p7_adapter_extension), continuing on past the raw (P7-free) DNA+P5 column -- drawn at the
      same length as in draw_dna_fragments (just not necessarily the same x position, since
      excluded fragments shift everything after them) -- with the P7 pieces extending further
      out the bottom and/or top in this one.
    """
    pos_x, pos_y = position
    row_left_x = pos_x - dna_fragments_with_p7_row_width() / 2

    column_index = 0
    for i, (top_spans, bottom_spans, label_left, label_right) in enumerate(build_dna_fragments()):
        if i in EXCLUDED_FRAGMENT_INDICES:
            continue

        column_x = row_left_x + column_index * FRAGMENT_ROW_SPACING
        column_index += 1
        right_x = column_x + STRAND_SEPARATION / 2 + STRAND_WIDTH / 2
        left_x = column_x - STRAND_SEPARATION / 2 - STRAND_WIDTH / 2
        fragment_length = sum(length for length, _, _ in top_spans)
        bottom_y = pos_y - fragment_length / 2
        top_y = pos_y + fragment_length / 2

        draw_dna_fragment_column(left_x, bottom_y, top_spans)
        draw_dna_fragment_column(right_x, bottom_y, bottom_spans)

        # Antiparallel strands: a plain bottom end has the left (top-strand) column at 5', a
        # plain top end has the right (bottom-strand) column at 5' -- that's which column each
        # P7 adapter's orange tail sits on (mirrored from draw_dna_fragments' column assignment,
        # see there).
        if label_left:
            for is_top, x in ((True, left_x), (False, right_x)):
                p7_spans = p7_adapter_extension(is_tail=is_top, reversed_order=True)
                p7_length = sum(length for length, _, _ in p7_spans)
                draw_dna_fragment_column(x, bottom_y - p7_length, p7_spans)
        if label_right:
            for is_top, x in ((True, left_x), (False, right_x)):
                draw_dna_fragment_column(x, top_y, p7_adapter_extension(is_tail=not is_top, reversed_order=False))


def first_p7_fragment_tail_adapter_center():
    """
    Returns (x, y) of the center of the tail (dark magenta + dark orange) P7 adapter on the
    first fragment draw_dna_fragments_with_p7_adapters draws, relative to its own position --
    i.e. the same left_x/bottom_y geometry that function uses internally for its first
    (non-excluded) column, factored out so __main__ can point a label at it without duplicating
    that layout logic.
    """
    row_left_x = -dna_fragments_with_p7_row_width() / 2
    top_spans, _, _, _ = next(
        fragment for i, fragment in enumerate(build_dna_fragments()) if i not in EXCLUDED_FRAGMENT_INDICES
    )
    left_x = row_left_x - STRAND_SEPARATION / 2 - STRAND_WIDTH / 2
    fragment_length = sum(length for length, _, _ in top_spans)
    bottom_y = -fragment_length / 2
    p7_length = sum(length for length, _, _ in p7_adapter_extension(is_tail=True, reversed_order=True))
    return left_x, bottom_y - p7_length / 2


FLOW_CELL_LINE_LENGTH = 2.4
FLOW_CELL_LINE_THICKNESS = 0.06
FLOW_CELL_NUM_PRIMERS = 8
LIGATED_PRIMER_INDICES = (0, 2)
# Which primer the denatured strands' group starts to the right of -- kept as its own constant,
# independent of LIGATED_PRIMER_INDICES, so the two groups can be positioned (and shifted)
# separately instead of the denatured group always trailing directly behind the ligated one.
OTHER_STRANDS_ANCHOR_PRIMER_INDEX = 6
DENATURED_STRAND_SPACING = 0.15
# How far the denatured strands' group -- the 6 unattached vertical DNA strands drawn beside the
# flow cell -- sits to the right of the flow cell's OTHER_STRANDS_ANCHOR_PRIMER_INDEX primer, and
# above the primer tops -- set independently since the two aren't related (e.g. shrinking the
# vertical gap shouldn't pull the group leftward). This is the single horizontal-position knob
# for that whole group: the washed_away_label ("fragments that cannot bind...") in draw_flow_cell
# is placed relative to the strands' own drawn position, not given a separate offset, so it moves
# right along with them whenever this changes.



def find_span_center_offset(spans, color):
    """
    Returns the x offset (from the start of spans, i.e. the strand's own left edge) of the
    center of the first span whose color is `color`, or None if no span has that color.
    """
    offset = 0
    for length, facecolor, _ in spans:
        if facecolor == color:
            return offset + length / 2
        offset += length
    return None


def find_span_length(spans, color):
    """Returns the length of the (first) span in spans whose color is `color`."""
    for length, facecolor, _ in spans:
        if facecolor == color:
            return length
    return None


def draw_dna_fragment_column(x, y_start, spans):
    """
    Draws one column (one strand) of a standalone, closed DNA fragment, starting at y_start and
    stacking upward, from an ordered list of (length, facecolor, seamless_before) segments
    (bottom to top). Every segment is fully bordered, except where a segment's seamless_before
    is True, which leaves both sides of that seam borderless (the same physical piece) -- the
    vertical mirror of draw_dna_fragment_row.
    """
    y = y_start
    for i, (length, facecolor, seamless_before) in enumerate(spans):
        skip_edges = []
        if seamless_before:
            skip_edges.append("bottom")
        if i + 1 < len(spans) and spans[i + 1][2]:
            skip_edges.append("top")
        draw_strand_segment(x - STRAND_WIDTH / 2, y, STRAND_WIDTH, length, skip_edges=tuple(skip_edges), facecolor=facecolor)
        y += length


def split_ligated_and_other_strands():
    """
    Splits every strand from every (non-excluded) fragment build_dna_fragments builds -- each
    with a P7 adapter added on whichever end doesn't already have a P5 adapter, same as
    draw_dna_fragments_with_p7_adapters (see add_p7_adapters) -- into the (at most two) that
    have a light green section and every other strand (denatured, unattached).

    The two with a light green section get ligated to the flow cell; since that end already has
    its P5 adapter, their new P7 adapter necessarily falls on their *other* (previously plain)
    end. That end's DNA-and-P5 spans aren't relevant to where the P7 part is drawn, so for these
    two, only the P7 extension itself is returned (in near-DNA-to-far order, ready to continue
    past their DNA in draw_flow_cell's vertical layout) alongside their original (P7-free) P5
    spans -- as (p5_spans, p7_extension) pairs. Every other strand is returned as a single flat
    spans list (P7 already merged in), ready to draw as-is. Returns (ligated_strands,
    other_strands).
    """
    ligated_strands = []
    other_strands = []
    for i, (top_spans, bottom_spans, label_left, label_right) in enumerate(build_dna_fragments()):
        if i in EXCLUDED_FRAGMENT_INDICES:
            continue
        for is_top, spans in ((True, top_spans), (False, bottom_spans)):
            if find_span_center_offset(spans, ADAPTER_COLOR_GREEN_LIGHT) is not None and len(ligated_strands) < 2:
                is_tail = is_top if label_left else not is_top
                ligated_strands.append((spans, p7_adapter_extension(is_tail)))
            else:
                other_strands.append(add_p7_adapters(spans, is_top, label_left, label_right))
    return ligated_strands, other_strands


def flow_cell_primer_tops():
    """Returns a dict of primer index -> its top's height above the flow cell's own y position."""
    p5_adapter_length = (STRAND_OVERHANG + DSB_WIDTH + 2 * FRAGMENT_SHIFT - ADAPTER_SET_GAP) / 2
    green_length = p5_adapter_length / 2
    orange_length = green_length
    return {
        i: FLOW_CELL_LINE_THICKNESS / 2 + (green_length if i % 2 == 0 else orange_length)
        for i in range(FLOW_CELL_NUM_PRIMERS)
    }


def flow_cell_content_height():
    """
    Returns the total height of draw_flow_cell's content, from its own y position (the black
    line) up to the top of its highest strand -- so callers positioning multiple diagrams (e.g.
    __main__) can space around it correctly without having to duplicate its layout logic.
    """
    primer_tops = flow_cell_primer_tops()
    ligated_strands, other_strands = split_ligated_and_other_strands()
    ligated_tops = [
        FLOW_CELL_LINE_THICKNESS / 2 + find_span_length(spans, ADAPTER_COLOR_GREEN_LIGHT)
        + find_span_length(spans, ADAPTER_COLOR_BLUE_LIGHT) + find_span_length(spans, "lightgray")
        + sum(length for length, _, _ in p7_extension)
        for spans, p7_extension in ligated_strands
    ]
    # The "other" strands sit beside the ligated ones now (see draw_flow_cell), not stacked
    # above them, so their own height only needs to clear the plain primers, not the (much
    # taller) ligated strands.
    other_max_length = max(sum(length for length, _, _ in spans) for spans in other_strands)
    other_top = max(primer_tops.values()) + DENATURED_STRANDS_VERTICAL_OFFSET + other_max_length
    return max(list(primer_tops.values()) + ligated_tops + [other_top])


def flow_cell_content_right_extent():
    """
    Returns how far draw_flow_cell's content reaches to the right of its own x position -- i.e.
    past FLOW_CELL_LINE_LENGTH / 2, since the "other" (denatured) strands sit further right than
    the line's own right end, beside the ligated strands. Used to size __main__'s xlim so that
    group doesn't get clipped.
    """
    _, other_strands = split_ligated_and_other_strands()
    primer_spacing = FLOW_CELL_LINE_LENGTH / (FLOW_CELL_NUM_PRIMERS + 1)
    anchor_x = -FLOW_CELL_LINE_LENGTH / 2 + (OTHER_STRANDS_ANCHOR_PRIMER_INDEX + 1) * primer_spacing + STRAND_WIDTH
    other_x = anchor_x + STRAND_WIDTH / 2 + DENATURED_STRANDS_HORIZONTAL_OFFSET
    other_x += (len(other_strands) - 1) * DENATURED_STRAND_SPACING
    # + 1.1: the "washed away" label drawn just right of this group, wrapped to that same max
    # width (see draw_flow_cell) -- a safe upper bound on how far it reaches without having to
    # duplicate the label's own measured-text positioning logic here.
    return other_x + STRAND_WIDTH / 2 + ARROW_LABEL_OFFSET + 0.7


def draw_flow_cell(position=(0, 0)):
    """
    Draws a thick black horizontal line (a flow cell surface) centered at position=(x, y), with
    FLOW_CELL_NUM_PRIMERS primer rectangles evenly spaced along it and sticking straight up,
    perpendicular to the line. They alternate dark green and dark orange, all the same shape and
    size: that of the dark green section of a P5 adapter (see build_dna_fragments), rotated 90
    degrees to stick up instead of lying sideways along a DNA strand.

    Two individual (single, not paired) DNA strands are then shown ligated to two of the dark
    green primers: whichever strand, out of every fragment build_dna_fragments builds (skipping
    EXCLUDED_FRAGMENT_INDICES, matching draw_dna_fragments_with_p7_adapters), has a light green
    section -- there are always exactly two, one each from the two P5-adapter-bearing fragments,
    since only one of their two strands has the light (3') shade. Each is drawn as a vertical
    column standing beside (not on top of) one of the primers at LIGATED_PRIMER_INDICES: its
    light green section sits at the very bottom, at the same height as the primer and sharing
    its long (right) edge, then its light blue section continues seamlessly upward from the top
    of that, then its gray DNA section (a normal, bordered junction) extending furthest up --
    i.e. the two-tone adapter sits nearest the surface and the plain DNA points away from it,
    regardless of which order those sections happen to be built in horizontally elsewhere. Since
    that end already carries its P5 adapter, each of these strands' *other* end -- always its 5'
    terminus, since a strand has exactly one -- gets a P7 adapter too (see p7_adapter_extension),
    continuing on past the top of the gray DNA section.

    Every other strand from those fragments (the denatured, unattached ones) is drawn beside the
    flow cell too, each as its own vertical column, side by side with DENATURED_STRAND_SPACING
    between every one of them -- a single, uniform gap, unlike the paired diagrams elsewhere in
    this file where the within-pair spacing (STRAND_SEPARATION) differs from the
    between-fragment spacing (FRAGMENT_ROW_SPACING). To conserve space, this group sits beside
    the flow cell rather than stacked above it, starting to the right of the primer at
    OTHER_STRANDS_ANCHOR_PRIMER_INDEX (independent of where the ligated strands sit, at
    LIGATED_PRIMER_INDICES) and from a baseline level with the primer tops (not the much taller
    ligated strands), so their (differing) lengths extend upward by differing amounts from that
    lower baseline. Every one of them also
    gets a P7 adapter on whichever end(s) don't already have a P5 adapter, exactly as in
    draw_dna_fragments_with_p7_adapters (see add_p7_adapters) -- so, unlike draw_dna_fragments,
    no strand in this diagram is left plain.
    """
    pos_x, pos_y = position

    line_left_x = pos_x - FLOW_CELL_LINE_LENGTH / 2
    draw_strand_segment(line_left_x, pos_y - FLOW_CELL_LINE_THICKNESS / 2, FLOW_CELL_LINE_LENGTH, FLOW_CELL_LINE_THICKNESS, facecolor="black")

    primer_spacing = FLOW_CELL_LINE_LENGTH / (FLOW_CELL_NUM_PRIMERS + 1)
    primer_top_y = {i: pos_y + height for i, height in flow_cell_primer_tops().items()}
    for i, top_y in primer_top_y.items():
        primer_x = line_left_x + (i + 1) * primer_spacing
        is_green = i % 2 == 0
        color = ADAPTER_COLOR_GREEN_DARK if is_green else ADAPTER_COLOR_ORANGE_DARK
        draw_strand_segment(primer_x - STRAND_WIDTH / 2, pos_y + FLOW_CELL_LINE_THICKNESS / 2, STRAND_WIDTH, top_y - (pos_y + FLOW_CELL_LINE_THICKNESS / 2), facecolor=color)

    ligated_strands, other_strands = split_ligated_and_other_strands()

    primer_base_y = pos_y + FLOW_CELL_LINE_THICKNESS / 2
    ligated_tops = []
    first_strand_x = None
    first_strand_top = None
    for (spans, p7_extension), primer_index in zip(ligated_strands, LIGATED_PRIMER_INDICES):
        # The strand sits immediately to the right of its primer, so the light green section
        # shares the primer's long (right) edge instead of sitting on top of it.
        strand_x = line_left_x + (primer_index + 1) * primer_spacing + STRAND_WIDTH
        green_length = find_span_length(spans, ADAPTER_COLOR_GREEN_LIGHT)
        blue_length = find_span_length(spans, ADAPTER_COLOR_BLUE_LIGHT)
        gray_length = find_span_length(spans, "lightgray")

        draw_strand_segment(strand_x - STRAND_WIDTH / 2, primer_base_y, STRAND_WIDTH, green_length, skip_edges=("top",), facecolor=ADAPTER_COLOR_GREEN_LIGHT)
        upper_spans = [
            (blue_length, ADAPTER_COLOR_BLUE_LIGHT, True),
            (gray_length, "lightgray", False),
        ] + p7_extension
        draw_dna_fragment_column(strand_x, primer_base_y + green_length, upper_spans)
        strand_top = primer_base_y + green_length + sum(length for length, _, _ in upper_spans)
        ligated_tops.append(strand_top)
        if first_strand_x is None:
            first_strand_x = strand_x
            first_strand_top = strand_top

    # Explain, to the left of the first ligated strand, why only some fragments end up ligated
    # to the flow cell at all -- wrapped to a fixed data-coordinate width (measured against the
    # actual rendered text, see wrap_label_to_width) rather than a guessed character count, so
    # it can't run into the strand regardless of the label's length or the figure's DPI/font.
    first_strand_label = (
        "Only fragments from a DSB can bind to the flow cell. "
        "These are sequenced, and from this the location of the DSB can be determined."
    )
    wrapped_first_strand_label = wrap_label_to_width(first_strand_label, 1.1, fontsize=NON_ARROW_LABEL_FONT_SIZE)
    ax.text(
        first_strand_x - STRAND_WIDTH / 2 - ARROW_LABEL_OFFSET, (primer_base_y + first_strand_top) / 2,
        wrapped_first_strand_label, ha="right", va="center", fontsize=NON_ARROW_LABEL_FONT_SIZE,
    )

    other_base_y = max(primer_top_y.values()) + DENATURED_STRANDS_VERTICAL_OFFSET
    anchor_x = line_left_x + (OTHER_STRANDS_ANCHOR_PRIMER_INDEX + 1) * primer_spacing + STRAND_WIDTH
    other_x = anchor_x + STRAND_WIDTH / 2 + DENATURED_STRANDS_HORIZONTAL_OFFSET
    other_max_length = max(sum(length for length, _, _ in spans) for spans in other_strands)
    for spans in other_strands:
        draw_dna_fragment_column(other_x, other_base_y, spans)
        other_x += DENATURED_STRAND_SPACING

    # Explain, to the right of the last unligated strand, why they're drawn detached from the
    # flow cell line -- wrapped the same way as first_strand_label, for the same reason.
    washed_away_label = "Fragments that cannot bind to the flow cell are washed away."
    wrapped_washed_away_label = wrap_label_to_width(washed_away_label, 0.7, fontsize=NON_ARROW_LABEL_FONT_SIZE)
    ax.text(
        other_x - DENATURED_STRAND_SPACING + STRAND_WIDTH / 2 + ARROW_LABEL_OFFSET, other_base_y + other_max_length / 2,
        wrapped_washed_away_label, ha="left", va="center", fontsize=NON_ARROW_LABEL_FONT_SIZE,
    )


LEGEND_SWATCH_WIDTH = 0.35
LEGEND_SWATCH_HEIGHT = STRAND_WIDTH
LEGEND_TEXT_GAP = 0.08
# Vertical gap between each consecutive pair of the legend's 5 elements (the title, then its 4
# rows), in order -- so LEGEND_ROW_SPACINGS[0] is the gap between the title and the DNA strand
# row, LEGEND_ROW_SPACINGS[1] between the DNA strand and P5 adapter rows, and so on. 5 elements
# means exactly 4 gaps between them, so draw_legend expects exactly 4 entries here.
LEGEND_ROW_SPACINGS = [0.12, 0.32, 0.4, 0.4]
LEGEND_TITLE_FONT_SIZE = LABEL_FONT_SIZE
LEGEND_DARK_SWATCH_GAP = 0.03
LEGEND_DARK_COLORS = [ADAPTER_COLOR_BLUE_DARK, ADAPTER_COLOR_GREEN_DARK, ADAPTER_COLOR_MAGENTA_DARK, ADAPTER_COLOR_ORANGE_DARK]
# How far above, left of, and below the legend's own content (the title, swatches, and labels
# drawn by draw_legend) its box is drawn -- these three set the box's top, left, and bottom edges
# respectively. The box's right edge isn't tied to the content at all: it's just LEGEND_BOX_WIDTH
# to the right of the left edge these padding constants place, so LEGEND_BOX_WIDTH needs to be
# wide enough to cover LEGEND_BOX_PADDING_LEFT + LEGEND_SWATCH_WIDTH + LEGEND_TEXT_GAP + the
# widest label, or that label will run past the box's right edge.
LEGEND_BOX_PADDING_ABOVE = 0.1
LEGEND_BOX_PADDING_LEFT = 0.05
LEGEND_BOX_PADDING_BELOW = 0.16
# Total width of the box drawn around the legend, from the left edge LEGEND_BOX_PADDING_LEFT
# places to the box's right edge -- the one knob for how wide it is (see the padding constants
# above for why it needs to be big enough for the content it's covering).
LEGEND_BOX_WIDTH = 1.2
LEGEND_BOX_LINEWIDTH = 1


def draw_legend(position=(0, 0)):
    """
    Draws the diagrams' legend inside a box: a "Legend" title followed by four entries stacked
    below it (spaced according to LEGEND_ROW_SPACINGS), every entry a swatch with its label just
    to the right of it, vertically centered on the same row:
    - a plain gray swatch, the same facecolor used for DNA elsewhere in this file, labeled
      "DNA strand".
    - a light P5 adapter swatch -- the same light-blue/light-green two-tone strand ligated onto
      DNA ends in the other diagrams -- labeled to call out that its green half is the sequence
      that binds the flow cell.
    - a light P7 adapter swatch -- light magenta/light orange -- labeled the same way for its
      orange half. Note the other diagrams never actually draw a light-orange P7 section: only a
      *tail* P7 adapter (see p7_adapter_extension) has an orange part, and a tail is always dark.
      This swatch is legend-only, shown in light shades for symmetry with the P5 entry above it.
    - a row of LEGEND_DARK_COLORS (blue, green, magenta, orange -- the dark shade of each
      adapter color used elsewhere), packed into the same LEGEND_SWATCH_WIDTH the other entries'
      single swatch occupies, with one shared label explaining that dark shading marks a
      complementary sequence wherever it appears in the other diagrams.

    position=(x, y) anchors the content itself (not the box): x is the swatches' left edge, y is
    the "Legend" title's vertical center -- the same anchor draw_legend used before it grew a box
    around it. The box is then derived from that content, using the LEGEND_BOX_PADDING_* and
    LEGEND_BOX_WIDTH constants documented above.

    Returns the y position of the box's bottom edge, so __main__ can size the space it needs.
    """
    swatch_left_x, title_y = position
    label_x = swatch_left_x + LEGEND_SWATCH_WIDTH + LEGEND_TEXT_GAP
    box_left_x = swatch_left_x - LEGEND_BOX_PADDING_LEFT
    box_center_x = box_left_x + LEGEND_BOX_WIDTH / 2

    y = title_y
    ax.text(box_center_x, y, "Legend", ha="center", va="center", fontsize=LEGEND_TITLE_FONT_SIZE, fontweight="bold")
    y -= LEGEND_ROW_SPACINGS[0]

    draw_strand_segment(swatch_left_x, y - LEGEND_SWATCH_HEIGHT / 2, LEGEND_SWATCH_WIDTH, LEGEND_SWATCH_HEIGHT, facecolor="lightgray")
    ax.text(label_x, y, "DNA strand", ha="left", va="center", fontsize=NON_ARROW_LABEL_FONT_SIZE)
    y -= LEGEND_ROW_SPACINGS[1]

    draw_two_tone_strand(swatch_left_x, y - LEGEND_SWATCH_HEIGHT / 2, LEGEND_SWATCH_WIDTH, LEGEND_SWATCH_HEIGHT, ADAPTER_COLOR_BLUE_LIGHT, ADAPTER_COLOR_GREEN_LIGHT)
    ax.text(label_x, y, "P5 adapter;\nsequence for\nbinding to flow\ncell in green", ha="left", va="center", fontsize=NON_ARROW_LABEL_FONT_SIZE)
    y -= LEGEND_ROW_SPACINGS[2]

    draw_two_tone_strand(swatch_left_x, y - LEGEND_SWATCH_HEIGHT / 2, LEGEND_SWATCH_WIDTH, LEGEND_SWATCH_HEIGHT, ADAPTER_COLOR_MAGENTA_LIGHT, ADAPTER_COLOR_ORANGE_LIGHT)
    ax.text(label_x, y, "P7 adapter;\nsequence for\nbinding to flow\ncell in orange", ha="left", va="center", fontsize=NON_ARROW_LABEL_FONT_SIZE)
    y -= LEGEND_ROW_SPACINGS[3]

    dark_swatch_width = (LEGEND_SWATCH_WIDTH - (len(LEGEND_DARK_COLORS) - 1) * LEGEND_DARK_SWATCH_GAP) / len(LEGEND_DARK_COLORS)
    x = swatch_left_x
    for color in LEGEND_DARK_COLORS:
        draw_strand_segment(x, y - LEGEND_SWATCH_HEIGHT / 2, dark_swatch_width, LEGEND_SWATCH_HEIGHT, facecolor=color)
        x += dark_swatch_width + LEGEND_DARK_SWATCH_GAP
    ax.text(label_x, y, "dark colours\nrepresent\ncomplimentary\nsequences", ha="left", va="center", fontsize=NON_ARROW_LABEL_FONT_SIZE)

    box_top_y = title_y + LEGEND_BOX_PADDING_ABOVE
    box_bottom_y = y - LEGEND_SWATCH_HEIGHT / 2 - LEGEND_BOX_PADDING_BELOW
    draw_rectangle(box_left_x, box_bottom_y, LEGEND_BOX_WIDTH, box_top_y - box_bottom_y, facecolor="none", edgecolor="black", linewidth=LEGEND_BOX_LINEWIDTH)

    return box_bottom_y


def draw_dna():
    pass


def measure_text_width(text, fontsize=LABEL_FONT_SIZE):
    """
    Returns the actual rendered width (in data coordinates) of text at fontsize, drawn on the
    module-level axes. Font metrics vary by character, and data-to-pixel scale depends on the
    figure's DPI and axes limits, so this draws (and immediately discards) a probe text object
    rather than guessing from character count.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    probe = ax.text(0, 0, text, fontsize=fontsize)
    bbox = probe.get_window_extent(renderer=renderer)
    (x0, _), (x1, _) = ax.transData.inverted().transform([[bbox.x0, 0], [bbox.x1, 0]])
    probe.remove()
    return x1 - x0


def wrap_label_to_width(label, max_width, fontsize=LABEL_FONT_SIZE):
    """
    Wraps label onto as few lines as possible such that, once actually rendered at fontsize, no
    line is wider than max_width (in data coordinates). Unlike textwrap.wrap with a fixed
    characters-per-line guess, this measures the real rendered width of each candidate wrapping
    (see measure_text_width), so a label can't silently overflow into whatever's placed beside it.
    """
    for wrap_width in range(len(label), 0, -1):
        lines = textwrap.wrap(label, width=wrap_width)
        if max(measure_text_width(line, fontsize) for line in lines) <= max_width:
            return "\n".join(lines)
    return "\n".join(textwrap.wrap(label, width=1))


DIAGRAM_SPACING = 0.5
# Length of the arrow leading from diagram 4 into diagram 5, and how close its top (start) and
# bottom tip sit to diagram 4 and diagram 5 respectively -- sized to exactly fit
# DIAGRAM_4_TO_PAIR_GAP (below), with a larger top margin than bottom margin so the arrow's start
# sits noticeably below diagram 4 rather than hugging it.
DIAGRAM_4_TO_5_ARROW_LENGTH = 0.25
DIAGRAM_4_TO_5_ARROW_TOP_MARGIN = 0.15
DIAGRAM_4_TO_5_ARROW_BOTTOM_MARGIN = 0.02
DIAGRAM_4_TO_PAIR_GAP = DIAGRAM_4_TO_5_ARROW_TOP_MARGIN + DIAGRAM_4_TO_5_ARROW_LENGTH + DIAGRAM_4_TO_5_ARROW_BOTTOM_MARGIN
# Gap between the pair (diagrams 5/6) and the flow cell diagram, and how much empty margin is
# kept below the flow cell -- both kept separate from DIAGRAM_SPACING (and smaller than it) so
# everything below diagram 4 sits a bit more tightly packed, reducing the figure's overall height.
PAIR_TO_FLOW_CELL_GAP = 0.4
FLOW_CELL_BOTTOM_MARGIN = 0.3
# Space kept below the legend's last row, mirroring FLOW_CELL_BOTTOM_MARGIN's role above it.
LEGEND_BOTTOM_MARGIN = 0.2
# How far right of center diagram 7 (draw_flow_cell) is drawn -- the flow cell line itself, its
# ligated fragments, and both of its labels (first_strand_label and washed_away_label) are all
# positioned relative to draw_flow_cell's own x (its `position` argument), so shifting this one
# value moves all of them together.
FLOW_CELL_X_OFFSET = 0.1
# Horizontal gap between diagrams 5 and 6 -- kept separate from DIAGRAM_SPACING (and larger than
# it) so this one gap can be widened without affecting the vertical spacing between other
# diagrams, which also uses DIAGRAM_SPACING.
PAIR_HORIZONTAL_GAP = 0.9
# Position of the arrow from diagram 6 into diagram 7 (draw_flow_cell), and of its caption
# ("Binding of fragments to flow cell"), each as an (x, y) offset from their own default
# position -- the arrow's from (p7_fragments_x, the vertical center of the gap between the pair
# and the flow cell), the caption's from its usual spot just left of the arrow (see
# draw_labeled_arrow's label_dx/label_dy). Kept as two separate pairs so the caption can be
# repositioned without dragging the arrow along with it.
DIAGRAM_6_TO_7_ARROW_X_OFFSET = 0
DIAGRAM_6_TO_7_ARROW_Y_OFFSET = 0
DIAGRAM_6_TO_7_LABEL_X_OFFSET = 0
DIAGRAM_6_TO_7_LABEL_Y_OFFSET = 0


if __name__ == "__main__":
    # draw_dna_fragments and draw_dna_fragments_with_p7_adapters now lay their fragments out as
    # columns (see draw_dna_fragments), so what used to be their vertical stack height is now
    # their horizontal row width, and what used to be a single fragment's width (irrelevant to
    # spacing rows apart) is now each column's height -- relevant to spacing this whole pair
    # away from the diagrams above and below it.
    fragments_row_width = dna_fragments_row_width()
    p7_row_width = dna_fragments_with_p7_row_width()
    pair_top_offset, pair_bottom_offset = dna_fragments_pair_extents()
    fragments_position_y = -3 * DIAGRAM_SPACING - DIAGRAM_4_TO_PAIR_GAP - pair_top_offset + BELOW_DIAGRAM_4_VERTICAL_SHIFT

    # The pair sits side by side (same y) rather than stacked. draw_dna_fragments_with_p7_adapters
    # skips EXCLUDED_FRAGMENT_INDICES and closes up the gap, so its row is narrower than
    # draw_dna_fragments' -- that freed width is added to the gap between the two (on top of the
    # usual PAIR_HORIZONTAL_GAP margin) so the arrow between them can stretch to fill it, rather
    # than going to waste as uneven whitespace.
    row_width_gap = fragments_row_width - p7_row_width
    pair_row_gap = PAIR_HORIZONTAL_GAP + row_width_gap
    fragments_x = -(pair_row_gap / 2 + fragments_row_width / 2)
    p7_fragments_x = pair_row_gap / 2 + p7_row_width / 2

    draw_dsb((0, 0))
    draw_blunted_dsb((0, -DIAGRAM_SPACING))
    draw_blunted_dsb_with_adapters((0, -2 * DIAGRAM_SPACING))
    draw_fragment_dna((0, -3 * DIAGRAM_SPACING), break_positions=BREAK_POSITIONS)
    draw_dna_fragments((fragments_x, fragments_position_y))
    draw_dna_fragments_with_p7_adapters((p7_fragments_x, fragments_position_y))
    # The outermost column of each row sticks out past the row's nominal edge by
    # column_half_width (it's centered on its own x, not flush with the row edge), so nothing
    # placed between the two rows -- the arrow or its caption -- can go past this width on
    # either side of center without overlapping one of them.
    column_half_width = STRAND_SEPARATION / 2 + STRAND_WIDTH
    pair_row_gap_clearance = pair_row_gap - 2 * column_half_width
    arrow_length = pair_row_gap_clearance - 2 * ARROW_LABEL_OFFSET
    pair_arrow_y = fragments_position_y + PAIR_ARROW_VERTICAL_OFFSET
    draw_horizontal_labeled_arrow(0, pair_arrow_y, length=arrow_length)
    wrapped_pair_arrow_label = wrap_label_to_width(PAIR_ARROW_LABEL, pair_row_gap_clearance)
    pair_arrow_label_text = ax.text(0, pair_arrow_y - ARROW_LABEL_OFFSET, wrapped_pair_arrow_label, ha="center", va="top", fontsize=LABEL_FONT_SIZE)

    # Point out the tail P7 adapter on diagram 6's first fragment with its own label, placed in
    # the gap between diagrams 5 and 6 -- below pair_arrow_label's own text, but still within the
    # pair's already-reserved bottom extent (no extra room needed below it for the flow cell
    # diagram), since neither diagram actually has content down at this y level.
    adapter_offset_x, adapter_offset_y = first_p7_fragment_tail_adapter_center()
    adapter_x, adapter_y = p7_fragments_x + adapter_offset_x, fragments_position_y + adapter_offset_y
    wrapped_half_functional_label = wrap_label_to_width(HALF_FUNCTIONAL_P7_LABEL, 1.0, fontsize=NON_ARROW_LABEL_FONT_SIZE)
    half_functional_label_x = fragments_x / 2
    half_functional_label_text = ax.text(
        half_functional_label_x, adapter_y, wrapped_half_functional_label, ha="center", va="center", fontsize=NON_ARROW_LABEL_FONT_SIZE,
    )

    fig.canvas.draw()
    half_functional_label_bbox = half_functional_label_text.get_window_extent(renderer=fig.canvas.get_renderer())
    label_edge_x, label_center_y = ax.transData.inverted().transform(
        (half_functional_label_bbox.x1, (half_functional_label_bbox.y0 + half_functional_label_bbox.y1) / 2)
    )
    # Extended past the tail adapter's own center (adapter_x) to STRAND_SEPARATION + half a
    # strand width further right -- the left edge of the neighboring light-pink (non-tail) P7
    # column -- so the line's tip just touches that rectangle instead of stopping short of it.
    adapter_pointer_x = adapter_x + STRAND_SEPARATION + STRAND_WIDTH / 2
    ax.annotate(
        "", xy=(adapter_pointer_x, adapter_y), xytext=(label_edge_x, label_center_y),
        arrowprops=dict(arrowstyle="-", color="black", linewidth=1),
    )

    flow_cell_height = flow_cell_content_height()
    flow_cell_position_y = fragments_position_y - pair_bottom_offset - flow_cell_height - PAIR_TO_FLOW_CELL_GAP
    draw_flow_cell((FLOW_CELL_X_OFFSET, flow_cell_position_y))

    # The legend sits below everything else, horizontally centered like the diagrams above it --
    # its box's top edge FLOW_CELL_BOTTOM_MARGIN below the flow cell, and its box's left edge
    # LEGEND_BOX_WIDTH/2 left of center. draw_legend's position is the content's own anchor, not
    # the box's, so both are worked back from the box position via the padding constants it uses
    # to derive the box from that content.
    legend_box_top_y = flow_cell_position_y - FLOW_CELL_BOTTOM_MARGIN
    legend_box_left_x = -LEGEND_BOX_WIDTH / 2
    legend_content_x = legend_box_left_x + LEGEND_BOX_PADDING_LEFT
    legend_content_y = legend_box_top_y - LEGEND_BOX_PADDING_ABOVE
    legend_bottom_y = draw_legend((legend_content_x, legend_content_y))

    # Each diagram's own (top, bottom) y-extent, top-to-bottom, used to find the empty gap
    # between one diagram and the next so an arrow can be centered in it. The first four are
    # thin (just a couple of strand rows) relative to DIAGRAM_SPACING, so their own position is
    # a fine stand-in for both their top and bottom. draw_dna_fragments and
    # draw_dna_fragments_with_p7_adapters share one entry, since they're side by side at the
    # same y rather than being separate rows.
    diagram_extents = [
        (0, 0),
        (-DIAGRAM_SPACING, -DIAGRAM_SPACING),
        (-2 * DIAGRAM_SPACING, -2 * DIAGRAM_SPACING),
        (-3 * DIAGRAM_SPACING, -3 * DIAGRAM_SPACING),
        (fragments_position_y + pair_top_offset, fragments_position_y - pair_bottom_offset),
        (flow_cell_position_y + flow_cell_height, flow_cell_position_y),
    ]
    for i in range(len(diagram_extents) - 1):
        gap_top = diagram_extents[i][1]
        gap_bottom = diagram_extents[i + 1][0]
        if i == 3:
            # This gap leads into the side-by-side pair rather than a single centered diagram,
            # so point the arrow straight down at diagram 5 specifically (on the left) instead
            # of down the middle, where it would land between the two. Biased toward the bottom
            # of the gap (rather than centered) so it clears diagram 4 above it. gap_bottom
            # carries BELOW_DIAGRAM_4_VERTICAL_SHIFT (it's derived from fragments_position_y),
            # but this arrow sits right below diagram 4 itself rather than belonging to the
            # shifted block, so that shift is subtracted back out -- the arrow stays put right
            # below diagram 4 while diagram 5 (and everything below it) moves.
            arrow_y = (gap_bottom - BELOW_DIAGRAM_4_VERTICAL_SHIFT) + DIAGRAM_4_TO_5_ARROW_BOTTOM_MARGIN + DIAGRAM_4_TO_5_ARROW_LENGTH / 2
            draw_labeled_arrow(fragments_x, arrow_y, ARROW_LABELS[i], length=DIAGRAM_4_TO_5_ARROW_LENGTH)
        elif i == 4:
            # This gap leads out of the side-by-side pair, so start the arrow at diagram 6
            # specifically (on the right) rather than down the middle -- with its label on the
            # left (facing back toward diagram 6) instead of the usual right side.
            draw_labeled_arrow(
                p7_fragments_x + DIAGRAM_6_TO_7_ARROW_X_OFFSET, (gap_top + gap_bottom) / 2 + DIAGRAM_6_TO_7_ARROW_Y_OFFSET,
                ARROW_LABELS[i], length=LAST_ARROW_LENGTH, label_side="left",
                label_dx=DIAGRAM_6_TO_7_LABEL_X_OFFSET, label_dy=DIAGRAM_6_TO_7_LABEL_Y_OFFSET,
            )
        elif i < 3:
            # Shifted left of center so their labels have more room before running into the
            # strand ends' 5'/3' labels off to the right.
            draw_labeled_arrow(-0.2, (gap_top + gap_bottom) / 2, ARROW_LABELS[i])
        else:
            draw_labeled_arrow(0, (gap_top + gap_bottom) / 2, ARROW_LABELS[i])

    # The outermost column's own half-width (it's centered on its x, so it sticks out a bit
    # past the row's nominal edge) needs to be included too, or its far edge gets clipped by
    # the axes' data limits.
    pair_left_edge = fragments_x - fragments_row_width / 2 - column_half_width
    pair_right_edge = p7_fragments_x + p7_row_width / 2 + column_half_width
    pair_half_width = max(-pair_left_edge, pair_right_edge)
    half_width = max(STRAND_LENGTH * 1.2, pair_half_width, FLOW_CELL_X_OFFSET + flow_cell_content_right_extent(), LEGEND_BOX_WIDTH / 2)
    ax.set_xlim(-half_width, half_width)
    ax.set_ylim(legend_bottom_y - LEGEND_BOTTOM_MARGIN, STRAND_LENGTH * 0.6)
    ax.set_aspect("equal")
    ax.axis("off")

    # Figure height grown to fit the legend added below the flow cell -- kept proportional to the
    # extra vertical data range (legend_position_y replaces what used to be the bottom of the
    # figure) so text-to-shape scale stays the same as before the legend was added, rather than
    # everything shrinking to squeeze the new content into the old height.
    fig.set_size_inches(10, 21 + 4 * (flow_cell_position_y - legend_bottom_y))
    fig.savefig(OUTPUT_JPEG_PATH, dpi=150, bbox_inches="tight", pad_inches=0.1)
    print(f"Saved diagram to {OUTPUT_JPEG_PATH}")

    # plt.show()