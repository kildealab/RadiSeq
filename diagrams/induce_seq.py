import os

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUTPUT_JPEG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "induce_seq_diagrams.jpg")

fig, ax = plt.subplots()

STRAND_WIDTH = 0.04
STRAND_SEPARATION = 0.04
STRAND_LENGTH = 1.
STRAND_OVERHANG = 0.4
DSB_WIDTH = 0.05
LABEL_OFFSET = 0.05
LABEL_VERTICAL_OFFSET = 0.05
BREAK_POSITIONS = [ 0.2, 0.95]
STRAND_EDGE_WIDTH = 1

OTHER_FRAGMENT_LENGTHS = [0.3, 0.6, 0.9]
FRAGMENT_ROW_SPACING = 0.25

# Labels for the arrows drawn between consecutive diagrams in __main__, top-to-bottom. Edit
# freely -- there's one arrow per gap between diagrams, so this needs exactly 5 entries (the
# draw_dna_fragments / draw_dna_fragments_with_p7_adapters pair counts as one, since they're
# drawn side by side rather than as separate rows -- see PAIR_ARROW_LABEL for the sideways arrow
# between those two specifically).
ARROW_LABELS = [
    "dsb end blunting",
    "P5 adapter ligation",
    "fragmentation",
    "",
    "test 6",
]

# Label for the sideways arrow between draw_dna_fragments and draw_dna_fragments_with_p7_adapters.
PAIR_ARROW_LABEL = "P7 adapter ligation and size filtering"


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
    ax.text(x + offset, y, label, ha=ha, va="center", fontsize=10)


ARROW_LENGTH = 0.2
ARROW_LABEL_OFFSET = 0.1


def draw_labeled_arrow(x, y_center, label, length=ARROW_LENGTH):
    """
    Draws a short downward-pointing arrow centered at (x, y_center), with label as text just to
    its right, vertically centered on the arrow. Used to mark the transition between
    consecutive diagrams in __main__.
    """
    ax.annotate(
        "", xy=(x, y_center - length / 2), xytext=(x, y_center + length / 2),
        arrowprops=dict(arrowstyle="->", color="black", linewidth=1.5),
    )
    ax.text(x + ARROW_LABEL_OFFSET, y_center, label, ha="left", va="center", fontsize=10)


def draw_horizontal_labeled_arrow(x_center, y, label, length=ARROW_LENGTH):
    """
    Draws a short rightward-pointing arrow centered at (x_center, y), with label as text just
    above it, horizontally centered on the arrow. The horizontal counterpart to
    draw_labeled_arrow, used to mark the transition between two side-by-side diagrams.
    """
    ax.annotate(
        "", xy=(x_center + length / 2, y), xytext=(x_center - length / 2, y),
        arrowprops=dict(arrowstyle="->", color="black", linewidth=1.5),
    )
    ax.text(x_center, y + ARROW_LABEL_OFFSET, label, ha="center", va="bottom", fontsize=10)


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


def draw_dna_fragments(position=(0, 0)):
    """
    Draws a stack of standalone, complete DNA fragments, each row above the previous one
    (spaced FRAGMENT_ROW_SPACING apart). Unlike the other diagrams in this file, every fragment
    here is fully bordered on all sides: each one is a separate, complete molecule, rather than
    a piece of a longer sequence continuing off frame.

    Two of the fragments have the same lengths and orientations as the two DNA pieces (each
    with one ligated P5 adapter) from draw_fragment_dna / draw_blunted_dsb_with_adapters -- just
    not at the same positions, since here they're independent rows in this stack. The remaining
    fragments are plain (no adapter), with lengths taken from the OTHER_FRAGMENT_LENGTHS global.
    """
    pos_x, pos_y = position

    for i, (top_spans, bottom_spans, label_left, label_right) in enumerate(build_dna_fragments()):
        row_y = pos_y - i * FRAGMENT_ROW_SPACING
        bottom_y = row_y - STRAND_SEPARATION / 2 - STRAND_WIDTH
        top_y = row_y + STRAND_SEPARATION / 2
        fragment_length = sum(length for length, _, _ in top_spans)
        left_x = pos_x - fragment_length / 2

        draw_dna_fragment_row(top_y, left_x, top_spans)
        draw_dna_fragment_row(bottom_y, left_x, bottom_spans)

        # Antiparallel strands: the top strand runs 5'->3' left to right, the bottom 3'->5'.
        # (Not labeled here -- only draw_dsb, the first diagram, gets 5'/3' text labels.)


P7_ADAPTER_LENGTH = 0.3
EXCLUDED_FRAGMENT_INDICES = [2]


def draw_p7_adapter(x, top_y, bottom_y, side, top_is_5prime):
    """
    Draws a P7 adapter ligated to a DNA end at x, extending further left (side='left') or right
    (side='right') from x. Both strands' pink (magenta) sections are the same length
    (P7_ADAPTER_LENGTH / 2): whichever strand is attached to that end's 5' terminus (top_y if
    top_is_5prime, else bottom_y) continues on past its pink half with an equal-length orange
    half (no border on the seam between them, the same physical piece), reaching a total length
    of P7_ADAPTER_LENGTH; the other strand has no orange part, so it's just the pink half on its
    own -- shorter than the 5'-terminus strand, ending where the orange part would have started.
    The 5'-terminus strand uses the dark magenta/orange shades, and the 3'-terminus (plain)
    strand uses the light magenta (pink) shade, so the two opposite strands read as obviously
    different.
    """
    half_length = P7_ADAPTER_LENGTH / 2
    plain_start = x if side == "right" else x - half_length
    tail_start = x if side == "right" else x - P7_ADAPTER_LENGTH

    tail_y = top_y if top_is_5prime else bottom_y
    plain_y = bottom_y if top_is_5prime else top_y

    draw_strand_segment(plain_start, plain_y, half_length, STRAND_WIDTH, facecolor=ADAPTER_COLOR_MAGENTA_LIGHT)
    if side == "right":
        draw_two_tone_strand(tail_start, tail_y, P7_ADAPTER_LENGTH, STRAND_WIDTH, ADAPTER_COLOR_MAGENTA_DARK, ADAPTER_COLOR_ORANGE_DARK)
    else:
        draw_two_tone_strand(tail_start, tail_y, P7_ADAPTER_LENGTH, STRAND_WIDTH, ADAPTER_COLOR_ORANGE_DARK, ADAPTER_COLOR_MAGENTA_DARK)


def p7_adapter_extension(is_tail, reversed_order=False):
    """
    Returns the (length, facecolor, seamless_before) segments of a P7 adapter attached to one
    end of a strand, in near-DNA-to-far order: just a light magenta half if not is_tail, or a
    dark magenta half followed (seamlessly) by a dark orange half if is_tail -- matching
    draw_p7_adapter's coloring. If reversed_order, returns them far-to-near instead (with the
    seamless_before flags recomputed to match), for prepending before a DNA span that extends
    to its right rather than appending after one that extends to its left.
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
    and/or label_right), matching draw_p7_adapter's antiparallel-strand rule for which end is
    the 5' (tail) one: at a left end the top strand is 5', at a right end the bottom strand is.
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


def draw_dna_fragments_with_p7_adapters(position=(0, 0)):
    """
    Draws the same fragment stack as draw_dna_fragments, except:
    - any fragment index in EXCLUDED_FRAGMENT_INDICES is skipped entirely; the remaining
      fragments keep the exact row positions (row_y = pos_y - i * FRAGMENT_ROW_SPACING, using
      their original index i) they have in draw_dna_fragments, rather than closing the gap.
    - every DNA end that doesn't already have a P5 adapter (i.e. every end draw_dna_fragments
      would otherwise give a plain 5'/3' text label) is instead ligated with a P7 adapter (see
      draw_p7_adapter): two same-length strands, one of which is two-toned pink/orange.
    """
    pos_x, pos_y = position

    for i, (top_spans, bottom_spans, label_left, label_right) in enumerate(build_dna_fragments()):
        if i in EXCLUDED_FRAGMENT_INDICES:
            continue

        row_y = pos_y - i * FRAGMENT_ROW_SPACING
        bottom_y = row_y - STRAND_SEPARATION / 2 - STRAND_WIDTH
        top_y = row_y + STRAND_SEPARATION / 2
        fragment_length = sum(length for length, _, _ in top_spans)
        left_x = pos_x - fragment_length / 2
        right_x = left_x + fragment_length

        draw_dna_fragment_row(top_y, left_x, top_spans)
        draw_dna_fragment_row(bottom_y, left_x, bottom_spans)

        # Antiparallel strands: a plain left end has its top strand at 5', a plain right end
        # has its bottom strand at 5' -- that's the strand each P7 adapter's orange tail sits on.
        if label_left:
            draw_p7_adapter(left_x, top_y, bottom_y, side="left", top_is_5prime=True)
        if label_right:
            draw_p7_adapter(right_x, top_y, bottom_y, side="right", top_is_5prime=False)


FLOW_CELL_LINE_LENGTH = 2.4
FLOW_CELL_LINE_THICKNESS = 0.06
FLOW_CELL_NUM_PRIMERS = 8
LIGATED_PRIMER_INDICES = (2, 4)
DENATURED_STRAND_SPACING = FRAGMENT_ROW_SPACING


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
    return max(list(primer_tops.values()) + ligated_tops) + len(other_strands) * DENATURED_STRAND_SPACING


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

    Every other strand from those fragments (the denatured, unattached ones) is drawn above the
    flow cell too, each as its own single row, stacked with DENATURED_STRAND_SPACING between
    every one of them -- a single, uniform gap, unlike the paired diagrams elsewhere in this
    file where the within-pair spacing (STRAND_SEPARATION) differs from the between-fragment
    spacing (FRAGMENT_ROW_SPACING). Every one of them also gets a P7 adapter on whichever end(s)
    don't already have a P5 adapter, exactly as in draw_dna_fragments_with_p7_adapters (see
    add_p7_adapters) -- so, unlike draw_dna_fragments, no strand in this diagram is left plain.
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
        ligated_tops.append(primer_base_y + green_length + sum(length for length, _, _ in upper_spans))

    other_y = max(list(primer_top_y.values()) + ligated_tops) + DENATURED_STRAND_SPACING
    for spans in other_strands:
        fragment_length = sum(length for length, _, _ in spans)
        draw_dna_fragment_row(other_y, pos_x - fragment_length / 2, spans)
        other_y += DENATURED_STRAND_SPACING


def draw_dna():
    pass

DIAGRAM_SPACING = 0.5

if __name__ == "__main__":
    fragments_position_y = -STRAND_LENGTH * 0.6 - 3 * DIAGRAM_SPACING - DIAGRAM_SPACING
    num_fragments = len(OTHER_FRAGMENT_LENGTHS) + 2
    fragments_stack_height = (num_fragments - 1) * FRAGMENT_ROW_SPACING

    # draw_dna_fragments and draw_dna_fragments_with_p7_adapters sit side by side (same y) rather
    # than stacked, spaced far enough apart (half of each one's own widest fragment, plus a gap)
    # that neither's content can cross into the other's.
    pair_center_gap = dna_fragments_max_width() / 2 + dna_fragments_with_p7_max_width() / 2 + DIAGRAM_SPACING
    fragments_x = -pair_center_gap / 2
    p7_fragments_x = pair_center_gap / 2

    draw_dsb((0, 0))
    draw_blunted_dsb((0, -DIAGRAM_SPACING))
    draw_blunted_dsb_with_adapters((0, -2 * DIAGRAM_SPACING))
    draw_fragment_dna((0, -3 * DIAGRAM_SPACING), break_positions=BREAK_POSITIONS)
    draw_dna_fragments((fragments_x, fragments_position_y))
    draw_dna_fragments_with_p7_adapters((p7_fragments_x, fragments_position_y))
    draw_horizontal_labeled_arrow(0, fragments_position_y + FRAGMENT_ROW_SPACING / 2, PAIR_ARROW_LABEL)

    flow_cell_height = flow_cell_content_height()
    flow_cell_position_y = fragments_position_y - fragments_stack_height - flow_cell_height - DIAGRAM_SPACING
    draw_flow_cell((0, flow_cell_position_y))

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
        (fragments_position_y, fragments_position_y - fragments_stack_height),
        (flow_cell_position_y + flow_cell_height, flow_cell_position_y),
    ]
    for i in range(len(diagram_extents) - 1):
        gap_top = diagram_extents[i][1]
        gap_bottom = diagram_extents[i + 1][0]
        draw_labeled_arrow(0, (gap_top + gap_bottom) / 2, ARROW_LABELS[i])

    pair_half_width = pair_center_gap / 2 + max(dna_fragments_max_width(), dna_fragments_with_p7_max_width()) / 2
    ax.set_xlim(-max(STRAND_LENGTH * 1.2, pair_half_width), max(STRAND_LENGTH * 1.2, pair_half_width))
    ax.set_ylim(flow_cell_position_y - DIAGRAM_SPACING, STRAND_LENGTH * 0.6)
    ax.set_aspect("equal")

    fig.set_size_inches(8, 21)
    fig.savefig(OUTPUT_JPEG_PATH, dpi=150)
    print(f"Saved diagram to {OUTPUT_JPEG_PATH}")

    plt.show()