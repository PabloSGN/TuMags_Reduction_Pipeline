import sys
from os import path, walk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox
from astropy.io import fits as pyfits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import RangeSlider

def view_data(file_ID, flip_order=False):

    stokes_limits = {
        0: (0.2, 1.2),     # Stokes I
        1: (-0.01, 0.01),  # Stokes Q
        2: (-0.01, 0.01),  # Stokes U
        3: (-0.01, 0.01)   # Stokes V
    }

    class CursorApp:
        def __init__(self, ax, stokes, axes, lines, wave_axis, im_ax1, clim=[0, 2]):
            self.ax = ax
            self.stokes = stokes
            self.ax2, self.ax3, self.ax4, self.ax5 = axes
            self.line2, self.line3, self.line4, self.line5 = lines
            self.wave_axis = wave_axis
            self.im_ax1 = im_ax1
            self.clim = clim

            self.lx = ax.axhline(color='red', linewidth=0.5)
            self.ly = ax.axvline(color='red', linewidth=0.5)
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.4)
            self.txt = ax.text(0, 0, '', fontsize=8, bbox=props)

            self.x = self.y = 0
            self.pressed = False
            self.stokes_index = 0
            self.wave_index = 0
            self.titles = ['I', 'Q', 'U', 'V']

        def mouse_move(self, event):
            if not event.inaxes:
                return
            if event.inaxes == self.ax and not self.pressed:
                self._mouse_move_event(event)

        def _mouse_move_event(self, event):
            sy, sx = self.stokes.shape[-2:]
            self.x = int(np.clip(round(event.xdata), 0, sx - 1))
            self.y = int(np.clip(round(event.ydata), 0, sy - 1))

            # Update text and crosshairs
            self.txt.set_position((self.x, self.y))
            self.txt.set_text(
                f'v = {self.stokes[self.wave_index, self.stokes_index,  self.y, self.x]:.4f}'
            )
            self.lx.set_ydata([self.y, self.y])
            self.ly.set_xdata([self.x, self.x])

            # Update right panel Stokes profiles
            self.line2.set_ydata(self.stokes[:, 0, self.y, self.x])
            self.line3.set_ydata(self.stokes[:, 1, self.y, self.x])
            self.line4.set_ydata(self.stokes[:, 2, self.y, self.x])
            self.line5.set_ydata(self.stokes[:, 3, self.y, self.x])

            # Update left panel dynamically
            self.im_ax1.set_data(self.stokes[self.wave_index, self.stokes_index, :, :])
            #self.im_ax1.set_clim(*self.clim)
            self.im_ax1.set_clim(*self.clim)

            # Redraw
            self.ax.figure.canvas.draw_idle()
            self.ax2.figure.canvas.draw_idle()
            self.ax3.figure.canvas.draw_idle()
            self.ax4.figure.canvas.draw_idle()
            self.ax5.figure.canvas.draw_idle()

        def on_press(self, event):
            if event.inaxes:
                self.pressed = True

        def on_release(self, event):
            self.pressed = False

        def key_press(self, event):
            if event.key == 'q':
                sys.exit(0)

        def on_click_zoom(self, event):
            """Permite hacer zoom con clic izquierdo (in) y derecho (out) en ax1."""
            if event.inaxes != self.ax:
                return

            # Factor de zoom (ajusta si quieres más/menos intensidad)
            zoom_factor = 1.5

            cur_xlim = self.ax.get_xlim()
            cur_ylim = self.ax.get_ylim()
            xdata, ydata = event.xdata, event.ydata

            if event.button == 1:  # Clic izquierdo = zoom in
                scale_factor = 1 / zoom_factor
            elif event.button == 3:  # Clic derecho = zoom out
                scale_factor = zoom_factor
            else:
                return  # Ignorar otros botones

            new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
            new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

            relx = (cur_xlim[1] - xdata) / (cur_xlim[1] - cur_xlim[0])
            rely = (cur_ylim[1] - ydata) / (cur_ylim[1] - cur_ylim[0])

            self.ax.set_xlim([xdata - new_width * (1 - relx),
                              xdata + new_width * relx])
            self.ax.set_ylim([ydata - new_height * (1 - rely),
                              ydata + new_height * rely])

            self.ax.figure.canvas.draw_idle()

    # -----------------------
    # Load FITS data
    # -----------------------
    with pyfits.open(file_ID, memmap=True) as hdu_list:
        stokes = hdu_list[0].data
    
    if stokes.ndim == 5:
        stokes = stokes[0]


    if flip_order:
        print("[INFO] Input order assumed to be (stokes, wave, y, x). Converting to (wave, stokes, y, x).")
        stokes = np.swapaxes(stokes, 0, 1)
    else:
        print("[INFO] Input order assumed to be (wave, stokes, y, x).")

    stokes = stokes / np.mean(stokes[-1, 0, :, :])
    wave_axis = np.arange(stokes.shape[0])

    # Si sólo hay I y V, crear Q y U artificiales
    if stokes.shape[1] == 2:

        stokes_iv = stokes

        stokes = np.zeros(
            (stokes.shape[0], 4, stokes.shape[2], stokes.shape[3]),
            dtype=stokes.dtype
        )

        stokes[:, 0] = stokes_iv[:, 0]  # I
        stokes[:, 1] = stokes_iv[:, 1]  # Q fake
        stokes[:, 2] = stokes_iv[:, 1]  # U fake
        stokes[:, 3] = stokes_iv[:, 1]  # V

        print("Detected NSTOKES=2 (I,V). Replicating V into Q,U,V for display.")
    # -----------------------
    # Create figure and axes
    # -----------------------
    fig = plt.figure(figsize=(14, 8))
    ax1 = fig.add_axes([0.05, 0.15, 0.4, 0.8])  # Left panel
    ax2 = fig.add_axes([0.55, 0.55, 0.18, 0.32])  # Q
    ax3 = fig.add_axes([0.75, 0.55, 0.18, 0.32])  # U
    ax4 = fig.add_axes([0.55, 0.20, 0.18, 0.32])  # V
    ax5 = fig.add_axes([0.75, 0.20, 0.18, 0.32])  # I

    # Sliders laterales para cada panel

    ax_rs2 = fig.add_axes([0.52, 0.55, 0.015, 0.32])
    ax_rs3 = fig.add_axes([0.72, 0.55, 0.015, 0.32])
    ax_rs4 = fig.add_axes([0.52, 0.20, 0.015, 0.32])
    ax_rs5 = fig.add_axes([0.72, 0.20, 0.015, 0.32])

    rs2 = RangeSlider(ax_rs2, "", 0.2, 1.2, valinit=(0.2, 1.2),
                    orientation="vertical")

    rs3 = RangeSlider(ax_rs3, "", -0.05, 0.05, valinit=(-0.01, 0.01),
                    orientation="vertical")

    rs4 = RangeSlider(ax_rs4, "", -0.05, 0.05, valinit=(-0.01, 0.01),
                    orientation="vertical")

    rs5 = RangeSlider(ax_rs5, "", -0.05, 0.05, valinit=(-0.01, 0.01),
                    orientation="vertical")

    def update_limits(val):

        ax2.set_ylim(rs2.val)
        ax3.set_ylim(rs3.val)
        ax4.set_ylim(rs4.val)
        ax5.set_ylim(rs5.val)

        fig.canvas.draw_idle()

    rs2.on_changed(update_limits)
    rs3.on_changed(update_limits)
    rs4.on_changed(update_limits)
    rs5.on_changed(update_limits)
    
    def force_right_axis(ax):
        # ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(True)
        ax.yaxis.set_ticks_position('right')
        ax.yaxis.set_label_position('right')
        ax.tick_params(left=False, labelleft=False)

    force_right_axis(ax3)
    force_right_axis(ax5)

    for ax in [ax2, ax3]:
        ax.tick_params(axis='x', labelbottom=False)

    fig.canvas.manager.set_window_title('Data inversion viewer')

    # -----------------------
    # Initialize left panel
    # -----------------------
    initial_clim = stokes_limits[0]
    im_ax1 = ax1.imshow(stokes[0, 0, :, :], cmap='inferno', clim=initial_clim, interpolation='none')

    # Save the original limits for reset
    original_xlim = ax1.get_xlim()
    original_ylim = ax1.get_ylim()

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im_ax1, cax=cax)
    ax1.set_title('Stokes 0, Wave 0')

    # -----------------------
    # Initialize right panel lines
    # -----------------------
    line2, = ax2.plot(wave_axis, stokes[:, 0, 0, 0], '.-')
    line3, = ax3.plot(wave_axis, stokes[:, 1, 0, 0], '.-')
    line4, = ax4.plot(wave_axis, stokes[:, 2, 0, 0], '.-')
    line5, = ax5.plot(wave_axis, stokes[:, 3, 0, 0], '.-')

    ax2.set_ylim((0.2, 1.2))
    ax3.set_ylim((-0.01, 0.01))
    ax4.set_ylim((-0.01, 0.01))
    ax5.set_ylim((-0.01, 0.01))

    # -----------------------
    # Bottom controls
    # -----------------------
    x0, x_gap, y_box, h_box, h_slider, y_slider = 0.05, 0.1, 0.04, 0.04, 0.04, 0.1
    x_right = 0.55

    # Stokes
    ax_box_stokes = fig.add_axes([x0, y_box, 0.05, h_box])
    box_stokes = TextBox(ax_box_stokes, 'Stokes', initial='0')
    ax_slider_stokes = fig.add_axes([x0, y_slider, 0.05, h_slider])
    slider_stokes = Slider(ax_slider_stokes, 'Stokes', 0, 3, valinit=0, valstep=1)
    box_stokes.label.set_fontsize(9)
    # box_stokes.label.set_horizontalalignment('right')
    slider_stokes.label.set_position((-0.15, 0.5))

    # Wavelength
    ax_box_wave = fig.add_axes([x0 + x_gap, y_box, 0.05, h_box])
    box_wave = TextBox(ax_box_wave, 'Wave', initial='0')
    ax_slider_wave = fig.add_axes([x0 + x_gap, y_slider, 0.05, h_slider])
    slider_wave = Slider(ax_slider_wave, 'Wave', 0, stokes.shape[0] - 1, valinit=0, valstep=1)
    box_wave.label.set_fontsize(9)
    # box_wave.label.set_horizontalalignment('right')
    slider_wave.label.set_position((-0.15, 0.5))

    # # Y-limits for right panel
    # ax_box_min = fig.add_axes([x_right, y_box+0.85, 0.05, h_box])
    # box_min = TextBox(ax_box_min, 'Y min', initial='-0.01')
    # ax_box_max = fig.add_axes([x_right + 0.2, y_box+0.85, 0.05, h_box])
    # box_max = TextBox(ax_box_max, 'Y max', initial='0.01')
    # box_min.label.set_fontsize(9)
    # # box_min.label.set_horizontalalignment('right')
    # box_max.label.set_fontsize(9)
    # # box_max.label.set_horizontalalignment('right')

    from matplotlib.widgets import Button

    # Reset Zoom button
    ax_button_reset = fig.add_axes([0.47, 0.05, 0.1, 0.05])  # [left, bottom, width, height]
    button_reset = Button(ax_button_reset, 'Reset Zoom', color='lightgray', hovercolor='0.9')

    # Auto Scale button
    ax_button_auto = fig.add_axes([0.59, 0.05, 0.1, 0.05])
    button_auto = Button(ax_button_auto, 'Auto Scale',
                        color='lightgray', hovercolor='0.9')

    def reset_zoom(event):
        ax1.set_xlim(original_xlim)
        ax1.set_ylim(original_ylim)
        fig.canvas.draw_idle()

    button_reset.on_clicked(reset_zoom)

    def auto_scale(event):

        x = cursor.x
        y = cursor.y

        # Escalas espectrales para el píxel seleccionado
        profI = stokes[:, 0, y, x]
        profQ = stokes[:, 1, y, x]
        profU = stokes[:, 2, y, x]
        profV = stokes[:, 3, y, x]


        def lims(x):
            p1, p99 = np.nanpercentile(x, [1, 99])
            if p1 == p99:
                p1 -= 1e-6
                p99 += 1e-6
            return p1, p99

        ax2.set_ylim(*lims(profI))
        ax3.set_ylim(*lims(profQ))
        ax4.set_ylim(*lims(profU))
        ax5.set_ylim(*lims(profV))

        # Reescala la imagen mostrada
        stokes_index = int(slider_stokes.val)
        wave_index = int(slider_wave.val)

        image = stokes[wave_index, stokes_index]

        im_ax1.set_clim(
            np.nanpercentile(image, 1),
            np.nanpercentile(image, 99)
        )

        cursor.clim = im_ax1.get_clim()

        fig.canvas.draw_idle()

    button_auto.on_clicked(auto_scale)
    # -----------------------
    # Cursor
    # -----------------------
    cursor = CursorApp(ax1, stokes, (ax2, ax3, ax4, ax5),
                       (line2, line3, line4, line5), wave_axis,
                       im_ax1,clim=initial_clim)
    fig.canvas.mpl_connect('motion_notify_event', cursor.mouse_move)
    fig.canvas.mpl_connect('key_press_event', cursor.key_press)
    fig.canvas.mpl_connect('button_press_event', cursor.on_press)
    fig.canvas.mpl_connect('button_release_event', cursor.on_release)
    fig.canvas.mpl_connect('button_press_event', cursor.on_click_zoom)

    # -----------------------
    # Callbacks
    # -----------------------
    def update_ax1(stokes_index=None, wave_index=None):
        if stokes_index is None:
            stokes_index = int(slider_stokes.val)
        if wave_index is None:
            wave_index = int(slider_wave.val)

        # Update CursorApp indices
        cursor.stokes_index = stokes_index
        cursor.wave_index = wave_index

        clim = stokes_limits[stokes_index]
        cursor.clim = clim  # synchronize with cursor too

        im_ax1.set_data(stokes[wave_index, stokes_index, :, :])
        im_ax1.set_clim(*clim)
        ax1.set_title(f'Stokes {stokes_index}, Wave {wave_index}')
        fig.canvas.draw_idle()

    def submit_stokes(text):
        try:
            val = int(text)
            val = np.clip(val, 0, 3)
            slider_stokes.set_val(val)
            update_ax1(stokes_index=val)
        except ValueError:
            pass

    def submit_wave(text):
        try:
            val = int(text)
            val = np.clip(val, 0, stokes.shape[0] - 1)
            slider_wave.set_val(val)
            update_ax1(wave_index=val)
        except ValueError:
            pass

    box_stokes.on_submit(submit_stokes)
    box_wave.on_submit(submit_wave)
    slider_stokes.on_changed(lambda val: update_ax1(stokes_index=int(val)))
    slider_wave.on_changed(lambda val: update_ax1(wave_index=int(val)))

    # # Y-limits
    # def submit_ymin(text):
    #     try:
    #         ymin = float(text)
    #         for ax in [ax3, ax4, ax5]:
    #             ax.set_ylim(bottom=ymin)
    #         # Update clim for the left panel if it matches current Stokes

    #         stokes_index = int(slider_stokes.val)
    #         old_clim = list(im_ax1.get_clim())
    #         im_ax1.set_clim(vmin=ymin, vmax=old_clim[1])
    #         cursor.clim = im_ax1.get_clim()  # sync with cursor
    #         fig.canvas.draw_idle()        
    #     except ValueError:
    #         pass

    # def submit_ymax(text):
    #     try:
    #         ymax = float(text)
    #         for ax in [ax3, ax4, ax5]:
    #             ax.set_ylim(top=ymax)
    #         # Update clim for the left panel if it matches current Stokes
    #         stokes_index = int(slider_stokes.val)
    #         old_clim = list(im_ax1.get_clim())
    #         im_ax1.set_clim(vmin=old_clim[0], vmax=ymax)
    #         cursor.clim = im_ax1.get_clim()
    #         fig.canvas.draw_idle()
    #     except ValueError:
    #         pass

    # box_min.on_submit(submit_ymin)
    # box_max.on_submit(submit_ymax)

    # -----------------------
    # Colormap selector
    # -----------------------
    from matplotlib.widgets import RadioButtons

    ax_color = fig.add_axes([0.30, 0.02, 0.1, 0.12])  # position on the figure
    cmap_options = ['inferno', 'plasma', 'viridis', 'gray', 'cividis', 'magma']
    radio_color = RadioButtons(ax_color, cmap_options, active=0)
    ax_color.set_title("Colormap", fontsize=9)

    def change_cmap(label):
        im_ax1.set_cmap(label)
        fig.canvas.draw_idle()

    radio_color.on_clicked(change_cmap)

    plt.show()


# ------------------------------------------------------------------
# Helper function unchanged
# ------------------------------------------------------------------
def list_fits(inpath: str = './', contain: str = None, remove_dir: bool = False,
              endswith=['.fits', '.fits.gz']):
    assert path.isdir(inpath)
    list_of_files = []
    for dirpath, _, filenames in walk(inpath):
        for filename in filenames:
            if filename.endswith(tuple(endswith)) and not filename.startswith('._'):
                filepath = path.join(dirpath, filename)
                if contain is None or contain in filename:
                    list_of_files.append(filename if remove_dir else filepath)
    return list_of_files


if __name__ == '__main__':

    def print_help():
        print("Usage: python visor.py <file.fits> [-flip]")
        print("Use -help to display this message.")
        print("")
        print("Options:")
        print("  -flip    Flip wavelength/order inside view_data")

    if len(sys.argv) >= 2 and sys.argv[1] == "-help":
        print_help()
        sys.exit(0)

    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Error: You must provide one FITS filename.")
        print_help()
        sys.exit(1)

    file_ID = sys.argv[1]

    flip_order = False

    if len(sys.argv) == 3:
        if sys.argv[2] == "-flip":
            flip_order = True
        else:
            print(f"Error: Unknown option {sys.argv[2]}")
            print_help()
            sys.exit(1)

    view_data(file_ID, flip_order=flip_order)