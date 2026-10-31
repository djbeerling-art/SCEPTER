# Print-size fonts: figures are placed 16.3 cm (6.4 in) wide in the Word file. On savefig, every text is scaled by
# K = figure width / 6.4 in (times fig.font_boost, default 1), so a font of N pt in the script prints at N pt;
# tick marks scale too, the last tight_layout call is redone, and the PNG is written at 300 dpi of print size.
import matplotlib.figure as mf, matplotlib.text as mt
PRINT_W = 6.4
BOOST = 0.9    # default: script font N pt prints at 0.9 N pt
_save, _tl = mf.Figure.savefig, mf.Figure.tight_layout
def _tight(self, *a, **k):
    self._tl_kw = k; return _tl(self, *a, **k)
def _savefig(self, *a, **k):
    if not getattr(self, '_scaled', False):
        K = self.get_figwidth()/PRINT_W*getattr(self, 'font_boost', BOOST)
        for t in self.findobj(mt.Text): t.set_fontsize(t.get_fontsize()*K)
        for ax in self.axes: ax.tick_params(length=3.5*K, width=0.8*K/1.3)
        if hasattr(self, '_tl_kw'): _tl(self, **self._tl_kw)
        self._scaled = True
    k['dpi'] = 300*PRINT_W/self.get_figwidth()
    return _save(self, *a, **k)
mf.Figure.savefig, mf.Figure.tight_layout = _savefig, _tight
