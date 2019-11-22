function images = import_images(exp,actions);
tic;
folder = exp.target_folder;
fprintf('Analyzing data from the %s folder,\n',folder)
import PD.Core.Extraction.get_load_tiffs;
[cutoutM] = get_load_tiffs(folder,actions.checkSameOrientation,actions.removeRegions,actions.makeSameSize);
 
import PD.Core.Extraction.register_images_fast;
images = register_images_fast(cutoutM);
t=toc;
fprintf('Image import completed in %.1f seconds. \n',t);
