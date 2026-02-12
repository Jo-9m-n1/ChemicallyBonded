# %%
import py3Dmol

view = py3Dmol.view(width=400, height=400, query='cid:222')

view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
view.zoomTo()
view.show()
