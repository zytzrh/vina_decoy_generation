import os
print(os.environ)

msms_path = "/NAS2020/Workspaces/DRLGroup/zyt/myenv/package/msms"
mgltools_path = "/NAS2020/Workspaces/DRLGroup/zyt/myenv/package/mgltools_x86_64Linux2_1.5.6"
vina_path = "/NAS2020/Workspaces/DRLGroup/zyt/myenv/package/autodock_vina_1_1_2_linux_x86/bin"
obabel_path = "/NAS2020/Workspaces/DRLGroup/zyt/myenv/Anaconda/anaconda3/envs/torchProtein38/bin"

# path for MSMS
os.environ["PATH"] += f":{msms_path}"
# set mgltool variable (if mac, should change mgltools_x86_64Linux2_1.5.6 into your downloaded mac version)
os.environ["PATH"] += f":{mgltools_path}/bin"
os.environ["MGL"] = f"{mgltools_path}/"
os.environ["MGLPY"] = f"{mgltools_path}/bin/python"
os.environ["MGLUTIL"] = f"{mgltools_path}/MGLToolsPckgs/AutoDockTools/Utilities24/"
# set vina dir  (if mac, should use /mac/release/)
os.environ["VINADIR"] = f"{vina_path}/"
# add path of obabel
os.environ["PATH"] += f":{obabel_path}/"

