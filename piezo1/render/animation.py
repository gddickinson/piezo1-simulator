"""Offscreen animation: render a sequence of frames and write a GIF or MP4.

An :class:`Animator` owns a standalone GL context and a scene, and steps a
user-supplied callback once per frame. The callback moves whatever it likes —
coordinates, camera, colours — and the animator captures the result.

Two choices here are about *legibility*, which is what separates a molecular
animation that teaches something from one that is merely busy:

* **Ease-in-out timing.** Linear motion between two states reads as a mechanical
  slide. A smoothstep makes the eye track the moving parts instead of the
  transition itself.
* **Held endpoints.** A held first and last frame gives the viewer time to see
  what the two states actually are, which is the whole point of showing a
  transition. Loops without it read as a jitter.

Captions are burned into the frame rather than left to the surrounding text,
because these files get pasted into slides and issue threads on their own.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

from ..config import RenderSettings
from .scene import Scene

__all__ = ["Animator", "AnimationSpec", "smoothstep", "ease_in_out",
           "ping_pong", "write_gif", "write_mp4"]


# --------------------------------------------------------------------------
# Timing
# --------------------------------------------------------------------------

def smoothstep(t: np.ndarray | float) -> np.ndarray | float:
    """Classic 3t² − 2t³ easing: zero velocity at both ends."""
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def ease_in_out(n_frames: int, hold: int = 0) -> np.ndarray:
    """Parameter values from 0 to 1 with eased timing and held endpoints."""
    moving = max(n_frames - 2 * hold, 2)
    t = smoothstep(np.linspace(0.0, 1.0, moving))
    return np.concatenate([np.zeros(hold), t, np.ones(hold)])


def ping_pong(n_frames: int, hold: int = 0) -> np.ndarray:
    """0 → 1 → 0, eased, with the endpoints held. Loops seamlessly."""
    half = ease_in_out(max(n_frames // 2, 2), hold=hold)
    return np.concatenate([half, half[::-1]])


# --------------------------------------------------------------------------
# Specification
# --------------------------------------------------------------------------

@dataclass
class AnimationSpec:
    """Everything needed to render one animation."""

    name: str
    n_frames: int = 60
    size: tuple[int, int] = (900, 700)
    fps: int = 20
    spin_degrees: float = 0.0
    caption: str = ""
    subcaption: str = ""
    samples: int = 4
    meta: dict = field(default_factory=dict)


# --------------------------------------------------------------------------
# Animator
# --------------------------------------------------------------------------

class Animator:
    """Renders frames offscreen into an in-memory list of images."""

    def __init__(self, spec: AnimationSpec,
                 settings: RenderSettings | None = None) -> None:
        import moderngl
        self.spec = spec
        self.ctx = moderngl.create_standalone_context(require=410)
        self.scene = Scene(self.ctx, settings or RenderSettings(samples=spec.samples))
        w, h = spec.size
        self.scene.resize(w, h)
        self._colour = self.ctx.texture((w, h), 4, samples=spec.samples)
        self._depth = self.ctx.depth_renderbuffer((w, h), samples=spec.samples)
        self._fbo = self.ctx.framebuffer(color_attachments=[self._colour],
                                         depth_attachment=self._depth)
        self._resolve_tex = self.ctx.texture((w, h), 4)
        self._resolve = self.ctx.framebuffer(color_attachments=[self._resolve_tex])
        self.frames: list = []

    # ------------------------------------------------------------- capture

    def capture(self):
        """Render the current scene state and return a PIL image."""
        from PIL import Image
        self._fbo.use()
        self.scene.render()
        self.ctx.finish()
        self.ctx.copy_framebuffer(self._resolve, self._fbo)
        w, h = self.spec.size
        data = self._resolve.read(components=3, alignment=1)
        return Image.frombytes("RGB", (w, h), data).transpose(Image.FLIP_TOP_BOTTOM)

    def run(self, step: Callable[[int, float], None],
            schedule: np.ndarray | None = None) -> list:
        """Render every frame, calling ``step(frame_index, t)`` beforehand."""
        sched = ease_in_out(self.spec.n_frames) if schedule is None else schedule
        self.frames = []
        for i, t in enumerate(sched):
            step(i, float(t))
            if self.spec.spin_degrees:
                self.scene.camera.spin(self.spec.spin_degrees / max(len(sched), 1))
            img = self.capture()
            if self.spec.caption:
                img = annotate(img, self.spec.caption, self.spec.subcaption,
                               progress=float(t))
            self.frames.append(img)
        return self.frames

    # -------------------------------------------------------------- output

    def save(self, path: str | Path, loop: bool = True,
             colors: int = 128, scale: float = 1.0) -> Path:
        """Write the captured frames.

        GIFs are downscaled and palette-reduced by default. A 960x720, 51-frame
        animation at full colour is ~8 MB, which is too heavy to live in a
        repository or paste into an issue; 0.75 scale with 128 colours brings
        that to roughly a quarter of the size with no visible loss on a dark
        background, where the palette is dominated by a few hues anyway.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix.lower() == ".gif":
            return write_gif(self.frames, path, fps=self.spec.fps, loop=loop,
                             colors=colors, scale=scale)
        return write_mp4(self.frames, path, fps=self.spec.fps)

    def release(self) -> None:
        for obj in (self._fbo, self._resolve, self._colour, self._depth,
                    self._resolve_tex):
            try:
                obj.release()
            except Exception:
                pass
        self.scene.release()
        try:
            self.ctx.release()
        except Exception:
            pass

    def __enter__(self) -> "Animator":
        return self

    def __exit__(self, *exc) -> None:
        self.release()


# --------------------------------------------------------------------------
# Frame annotation
# --------------------------------------------------------------------------

def annotate(img, caption: str, subcaption: str = "",
             progress: float | None = None):
    """Burn a caption and optional progress bar into a frame."""
    from PIL import ImageDraw
    d = ImageDraw.Draw(img)
    w, h = img.size
    d.text((16, 14), caption, fill=(232, 240, 255))
    if subcaption:
        d.text((16, 30), subcaption, fill=(150, 160, 178))
    if progress is not None:
        bar_w = w - 32
        y = h - 18
        d.rectangle([16, y, 16 + bar_w, y + 3], fill=(38, 44, 58))
        d.rectangle([16, y, 16 + bar_w * float(np.clip(progress, 0, 1)), y + 3],
                    fill=(91, 141, 239))
    return img


# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------

def write_gif(frames: list, path: str | Path, fps: int = 20,
              loop: bool = True, colors: int = 128,
              scale: float = 1.0) -> Path:
    """Write an animated GIF.

    GIF is limited to 256 colours, so an adaptive palette is generated from the
    whole sequence rather than per frame — a per-frame palette makes the
    background shimmer, which is far more distracting than the quantisation
    itself.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not frames:
        raise ValueError("no frames to write")
    from PIL import Image
    if scale != 1.0:
        w = int(frames[0].width * scale)
        h = int(frames[0].height * scale)
        frames = [f.resize((w, h), Image.LANCZOS) for f in frames]
    # Build one palette from a montage of sampled frames.
    sample = frames[:: max(len(frames) // 12, 1)]
    montage = Image.new("RGB", (frames[0].width, frames[0].height * len(sample)))
    for i, f in enumerate(sample):
        montage.paste(f, (0, i * frames[0].height))
    palette = montage.quantize(colors=colors, method=Image.MEDIANCUT)
    quantised = [f.quantize(palette=palette, dither=Image.FLOYDSTEINBERG)
                 for f in frames]
    quantised[0].save(path, save_all=True, append_images=quantised[1:],
                      duration=int(1000 / max(fps, 1)), loop=0 if loop else 1,
                      optimize=True)
    return path


def write_mp4(frames: list, path: str | Path, fps: int = 20,
              quality: int = 8) -> Path:
    """Write an MP4 via imageio/ffmpeg, falling back to GIF if unavailable."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import imageio.v2 as imageio
        writer = imageio.get_writer(str(path), fps=fps, quality=quality,
                                    macro_block_size=None)
        for f in frames:
            writer.append_data(np.asarray(f))
        writer.close()
        return path
    except Exception as exc:
        fallback = path.with_suffix(".gif")
        print(f"  MP4 unavailable ({type(exc).__name__}); writing {fallback.name}")
        return write_gif(frames, fallback, fps=fps)
