conda_env=mulled-v1-e05ad707e739a59dbca8c6a1fe3f0275ce0ad5649e8bab53f43146c97b48b37e
conda activate $conda_env
# Adjust threshold
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 255.0 'Default' 'red' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_default.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 255.0 'Percentile' 'over_under' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_percentile.gif' 'gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 255.0 'Huang' 'bw' 'yes' 'tools/image_processing/imagej2/test-data/blobs_threshold_huang_dark.gif' 'gif'
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 8.0 255.0 'Manual' 'bw' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_8-255.gif' 'gif'
# Test 5
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 8.0 'Manual' 'bw' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_0-8.gif' 'gif'
# Test 6
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel.tif' 0.0 255.0 'Default' 'bw' 'yes' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_threshold_default.tiff' 'tiff'


# Analyze particles
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Nothing' 'yes' 'no' 'no' '' '' 'tools/image_processing/imagej2/test-data/analyze_particles_nothing.tabular'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Outlines' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/analyze_particles_outlines.gif' 'gif' ''
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Bare_Outlines' 'no' 'yes' 'yes' 'tools/image_processing/imagej2/test-data/analyze_particles_bareoutlines.gif' 'gif' ''
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/analyze_particles_rois.tabular' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Overlay_Masks' 'yes' 'no' 'no' 'tools/image_processing/imagej2/test-data/analyze_particles_overlaymasks.gif' 'gif' 'tools/image_processing/imagej2/test-data/analyze_particles_nothing.tabular'
# Test 5
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/confocal_analyze_particles_rois.tabular' 'output_log.txt' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_threshold_default.tiff' 'yes' '0-Infinity' 0.0 1.0 'Overlay_Masks' 'yes' 'no' 'no' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_threshold_default_overlay_mask.tiff' 'tiff' 'tools/image_processing/imagej2/test-data/confocal_analyze_particles.tabular'

# Skeleton
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/skeletonized_blobs.gif' 'no' 'none' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/basic.tabular'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/skeletonized_clown.jpg' 'no' 'shortest_branch' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/shortest_branch_basic.tabular'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/skeletonized_blobs.gif' 'no' 'none' 'no' 'yes' 'no' '/tmp/tmpnbw7tzf2/job_working_directory/000/6/outputs/dataset_9ac5bcf4-4e41-4240-871e-7d110c52a7ff.dat'
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/skeletonized_clown.jpg' 'no' 'shortest_branch' 'no' 'yes' 'no' 'tools/image_processing/imagej2/test-data/shortest_branch_all_yes.tabular'

# binary_to_edm
# Test 1
ImageJ --ij2 --headless --debug --jython '/home/ldelisle/Documents/mygit/galaxytools/tools/image_processing/imagej2/imagej2_binary_to_edm_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 1 1 'no' 'no' 'tools/image_processing/imagej2/test-data/blobs_edm.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython '/home/ldelisle/Documents/mygit/galaxytools/tools/image_processing/imagej2/imagej2_binary_to_edm_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 10 3 'yes' 'yes' 'tools/image_processing/imagej2/test-data/blobs_black_edm.gif' 'gif'

# bunwarpj_adapt_transform
# Test 1
bunwarpj -adapt_transform 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/source_elastic_transformation.txt' 'tools/image_processing/imagej2/test-data/adapted_transformation.txt' 2.0

# bunwarpj_align
# Test 1
bunwarpj -align 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'NULL' 'tools/image_processing/imagej2/test-data/blobs.gif' 'NULL' 0 2 1 0.0 0.0 1.0 10.0 '/tmp/tmp8ch18obd/job_working_directory/000/7/outputs/dataset_520929f6-c649-427e-8606-15450fc1a39a.dat' 'tools/image_processing/imagej2/test-data/registered_source1.png' 'tools/image_processing/imagej2/test-data/registered_target1.png'
# Test 2
bunwarpj -align 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'NULL' 'tools/image_processing/imagej2/test-data/blobs.gif' 'NULL' 0 2 1 0.0 0.0 1.0 10.0 'tools/image_processing/imagej2/test-data/registered_source1.png' 'tools/image_processing/imagej2/test-data/registered_target1.png' '-save_transformation'
mv 'tools/image_processing/imagej2/test-data/registered_source1_transf.txt' 'tools/image_processing/imagej2/test-data/source_elastic_transformation_out_full.txt'
mv 'tools/image_processing/imagej2/test-data/registered_target1_transf.txt' 'tools/image_processing/imagej2/test-data/target_elastic_transformation_out_full.txt' 
# Test 3
bunwarpj -align 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/mask_white.png' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/mask_ramp.gif' 0 2 1 0.0 0.0 1.0 10.0 'tools/image_processing/imagej2/test-data/registered_source2.png' 'tools/image_processing/imagej2/test-data/registered_target2.png'

# bunwarpj_compare_elastic_raw
# Test 1
bunwarpj -compare_elastic_raw 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/target_elastic_transformation.txt' 'tools/image_processing/imagej2/test-data/source_raw_transformation.txt' 'tools/image_processing/imagej2/test-data/warping_index_raw_full.txt'  2>&1 | tee 'output_log.txt' && grep -Po 'Warping index = \K[^ ]+' 'output_log.txt' > 'tools/image_processing/imagej2/test-data/warping_index_raw_full.txt'

# bunwarpj....
# TODO

# imagej2...
# TODO

# Noise
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_noise_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'gif' 'add_specified_noise' 5.0 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/add_specified_noise.gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_noise_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'gif' 'despeckle' 'None' 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/despeckle.gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_noise_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'gif' 'remove_outliers' 'None' '10.0' '1' 'Bright' 'tools/image_processing/imagej2/test-data/remove_outliers.gif'
# Test 4 fails
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_noise_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'gif' 'remove_nans' 'None' 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/remove_nans.gif'

# Filter
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_filter_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'gif' 'gaussian_blur' 5.0 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/gaussian_blur.gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_filter_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'gif' 'median' 2.5 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/median.gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_filter_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'gif' 'unsharp_mask' 5.0 0.1 'None' 'None' 'tools/image_processing/imagej2/test-data/unsharp_mask.gif'
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_filter_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'gif' 'top_hat' 7.0 'None' 'light' 'dont' 'tools/image_processing/imagej2/test-data/top_hat.gif'
# Test 5
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_filter_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'gif' 'top_hat' 7.0 'None' '' '' 'tools/image_processing/imagej2/test-data/top_hat2.gif'
# Test 6
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_filter_jython_script.py' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel.tif'  'tiff' 'gaussian_blur' 2.0 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_gaussian_blur.tiff'

# Watershed
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_watershed_binary_jython_script.py' 'output.log' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' 'tools/image_processing/imagej2/test-data/blobs_watershed_binary.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_watershed_binary_jython_script.py' 'output.log' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_threshold_default.tiff' 'yes' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_default_threshold_watershed.tiff' 'tiff'
