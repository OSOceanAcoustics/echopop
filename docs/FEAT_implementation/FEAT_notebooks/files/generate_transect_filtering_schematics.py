import numpy as np
import matplotlib.pyplot as plt

# Color scheme
raw_color = "#888888"
raw_edge = "#222222"
on_color = "#3182bd"
on_edge = "#08519c"
off_color = "#e6550d"
off_edge = "#a63603"
filtered_color = "#31a354"
filtered_edge = "#006d2c"

# CASE 1: Trimming ends
# ---- Simulate points
nasc_pts = np.arange(0, 11, 1) # Raw NASC points/intervals from 0 to 10 nmi
on_effort = np.arange(1, 10, 1) # On-effort intervals from 1 to 9 nmi (inclusion filter)
# ---- Find overlapping points (apply inclusion filter)
off_effort_pts = nasc_pts[~np.isin(nasc_pts, on_effort)]
on_effort_pts = nasc_pts[np.isin(nasc_pts, on_effort)]
# ---- Plot
fig, ax = plt.subplots(figsize=(7, 2))

# Raw transect (gray, solid)
ax.plot(nasc_pts, np.full_like(nasc_pts, 1.1, dtype=float), color=raw_color, linewidth=5, zorder=1)
ax.scatter(nasc_pts, np.full_like(nasc_pts, 1.1, dtype=float), color=raw_color, edgecolor=raw_edge, s=90, marker="s", label="Raw transect", zorder=2)

# On-effort intervals (blue, solid)
ax.plot(on_effort, np.full_like(on_effort, 1.0), color=on_color, linewidth=5, zorder=3)
ax.scatter(on_effort, np.full_like(on_effort, 1.0), color=on_color, edgecolor=on_edge, s=90, marker="s", label="On-effort interval", zorder=4)
# Off-effort intervals (red, dashed) -- plot as segments
for pt in off_effort_pts:
    # For the left segment (0,1)
    if pt == off_effort_pts[0]:
        seg_x = [pt, pt+1]
    # For the right segment (9,10)
    elif pt == off_effort_pts[-1]:
        seg_x = [pt-1, pt]
    else:
        continue
    ax.plot(seg_x, [0.9, 0.9], color=off_color, linewidth=5, linestyle="--", zorder=5)
    ax.scatter(seg_x, [0.9, 0.9], color=off_color, edgecolor=off_edge, s=90, marker="s", label="Off-effort interval" if pt == off_effort_pts[0] else None, zorder=6)
    
    # Filtered output (green, dashed)
ax.plot(on_effort_pts, np.full_like(on_effort_pts, 0.8, dtype=float), color=filtered_color, linewidth=5, linestyle="-", zorder=7)
ax.scatter(on_effort_pts, np.full_like(on_effort_pts, 0.8, dtype=float), color=filtered_color, edgecolor=filtered_edge, s=90, marker="s", label="Filtered transect", zorder=8)
ax.set_ylim(0.7, 1.3)
ax.set_yticks([])
ax.set_xlabel("Distance (nmi)")
ax.legend(loc="upper center", ncol=4, fontsize=8.5)
plt.tight_layout()
plt.savefig("../../../_static/offeffort_transect_filtering_edges1.png", dpi=150, bbox_inches="tight")

plt.show()

# CASE 2: Multiple intervals trimmed at ends
# ---- Simulate points
on_effort = np.arange(2, 9, 1) # On-effort intervals from 2 to 8 nmi (inclusion filter)
# ---- Find overlapping points (apply inclusion filter)
off_effort_pts = nasc_pts[~np.isin(nasc_pts, on_effort)]
on_effort_pts = nasc_pts[np.isin(nasc_pts, on_effort)]
# ---- Plot
fig, ax = plt.subplots(figsize=(7, 2))

# Raw transect (gray, solid)
ax.plot(nasc_pts, np.full_like(nasc_pts, 1.1, dtype=float), color=raw_color, linewidth=5, zorder=1)
ax.scatter(nasc_pts, np.full_like(nasc_pts, 1.1, dtype=float), color=raw_color, 
           edgecolor=raw_edge, s=90, marker="s", label="Raw transect", zorder=2)

# On-effort intervals (blue, solid)
ax.plot(on_effort, np.full_like(on_effort, 1.0), color=on_color, linewidth=5, zorder=3)
ax.scatter(on_effort, np.full_like(on_effort, 1.0), color=on_color, edgecolor=on_edge, s=90, 
           marker="s", label="On-effort interval", zorder=4)
# Off-effort intervals (red, dashed) -- plot as segments
for pt in [0, 10]:
    # For the left segment (0,n)
    if pt == off_effort_pts[0]:
        seg_x = [pt, pt+1, pt+2]
    # For the right segment (10-n,10)
    elif pt == off_effort_pts[-1]:
        seg_x = [pt-2, pt-1, pt]
    else:
        continue
    ax.plot(seg_x, [0.9, 0.9, 0.9], color=off_color, linewidth=5, linestyle="--", zorder=5)
    ax.scatter(seg_x, [0.9, 0.9, 0.9], color=off_color, edgecolor=off_edge, s=90, marker="s", 
               label="Off-effort interval" if pt == off_effort_pts[0] else None, zorder=6)
    
    # Filtered output (green, dashed)
ax.plot(on_effort_pts, np.full_like(on_effort_pts, 0.8, dtype=float), color=filtered_color, 
        linewidth=5, linestyle="-", zorder=7)
ax.scatter(on_effort_pts, np.full_like(on_effort_pts, 0.8, dtype=float), color=filtered_color, 
           edgecolor=filtered_edge, s=90, marker="s", label="Filtered transect", zorder=8)
ax.set_ylim(0.7, 1.3)
ax.set_yticks([])
ax.set_xlabel("Distance (nmi)")
ax.legend(loc="upper center", ncol=4, fontsize=8.5)
plt.tight_layout()
plt.savefig("../../../_static/offeffort_transect_filtering_edges2.png", dpi=150, 
            bbox_inches="tight")

plt.show()

# CASE 3: Interior gaps
# ---- Simulate points
# -------- On-effort intervals with interior gaps
on_effort = np.concatenate([np.arange(0, 4), np.arange(5, 7), np.arange(8, 11)]) 
# ---- Find overlapping points (apply inclusion filter)
off_effort_pts = nasc_pts[~np.isin(nasc_pts, on_effort)]
on_effort_pts = nasc_pts[np.isin(nasc_pts, on_effort)]
# ---- Helper function to segment off-effort points
def get_segments(indices):
    """Return a list of contiguous segments from a sorted array of indices."""
    segments = []
    if len(indices) == 0:
        return segments
    start = indices[0]
    prev = indices[0]
    for idx in indices[1:]:
        if idx == prev + 1:
            prev = idx
        else:
            segments.append(np.arange(start, prev + 1))
            start = idx
            prev = idx
    segments.append(np.arange(start, prev + 1))
    return segments
# ---- Create segments
on_segments = get_segments(on_effort_pts)
off_segments = get_segments(off_effort_pts)
# ---- Plot
fig, ax = plt.subplots(figsize=(7, 2))

# Raw transect (gray, solid)
ax.plot(nasc_pts, np.full_like(nasc_pts, 1.1, dtype=float), color=raw_color, linewidth=5, zorder=1)
ax.scatter(nasc_pts, np.full_like(nasc_pts, 1.1, dtype=float), color=raw_color, edgecolor=raw_edge, 
           s=90, marker="s", label="Raw transect", zorder=2)

# On-effort intervals (blue, solid) -- segmented
for i, seg in enumerate(on_segments):
    ax.plot(seg, np.full_like(seg, 1.0, dtype=float), color=on_color, linewidth=5, zorder=3)
    ax.scatter(seg, np.full_like(seg, 1.0, dtype=float), color=on_color, edgecolor=on_edge, s=90, 
               marker="s", label="On-effort interval" if i == 0 else None, zorder=4)

# Off-effort intervals (orange-red, dashed) -- segmented
# Find the boundaries of each off-effort interval and plot a line between the adjacent on-effort 
# points
on_effort_sorted = np.sort(on_effort_pts)
gaps = np.where(np.diff(on_effort_sorted) > 1)[0]

for i, gap_idx in enumerate(gaps):
    start = on_effort_sorted[gap_idx]
    end = on_effort_sorted[gap_idx + 1]
    # Draw line between adjacent on-effort points
    ax.plot([start, end], [0.9, 0.9], color=off_color, linewidth=5, linestyle="--", zorder=5)
    # Scatter points for the off-effort interval (gap)
    off_pts = np.arange(start + 1, end)
    ax.scatter(off_pts, np.full_like(off_pts, 0.9, dtype=float), color=off_color, 
               edgecolor=off_edge, s=90, marker="s", 
               label="Off-effort interval" if i == 0 else None, zorder=6)
    # Scatter points for the ends of the interval
    ax.scatter([start, end], [0.9, 0.9], color=off_color, edgecolor=off_edge, s=90, marker="s", 
               zorder=6)
    
# Filtered output (green, dashed) -- segmented
for i, seg in enumerate(on_segments):
    ax.plot(seg, np.full_like(seg, 0.8, dtype=float), color=filtered_color, 
            linewidth=5, linestyle="-", zorder=7)
    ax.scatter(seg, np.full_like(seg, 0.8, dtype=float), color=filtered_color, 
               edgecolor=filtered_edge, s=90, marker="s", 
               label="Filtered transect" if i == 0 else None, zorder=8)

ax.set_ylim(0.7, 1.3)
ax.set_yticks([])
ax.set_xlabel("Distance (nmi)")
ax.legend(loc="upper center", ncol=4, fontsize=8.5)
plt.tight_layout()
plt.savefig("../../../_static/offeffort_transect_filtering_interior.png", dpi=150, 
            bbox_inches="tight")

plt.show()