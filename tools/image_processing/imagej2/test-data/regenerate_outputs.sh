source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda_env=mulled-v1-9e9e3735532d83dd44ff0f622055395263f964e81ca8113d4dfa7e133e3156bf
conda activate $conda_env
# Adjust threshold
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 255.0 'Default' 'red' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_default.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 255.0 'Percentile' 'over_under' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_percentile.gif' 'gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 255.0 'Huang' 'bw' 'yes' 'tools/image_processing/imagej2/test-data/blobs_threshold_huang_dark.gif' 'gif'
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 8.0 255.0 'Manual' 'bw' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_8-255.gif' 'gif'
# Test 5
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 8.0 'Manual' 'bw' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_0-8.gif' 'gif'
# Test 6
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel.tif' 0.0 255.0 'Default' 'bw' 'yes' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_threshold_default.tiff' 'tiff'


# Analyze particles
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Nothing' 'yes' 'no' 'no' '' '' 'tools/image_processing/imagej2/test-data/analyze_particles_nothing.tabular'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Outlines' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/analyze_particles_outlines.gif' 'gif' ''
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Bare_Outlines' 'no' 'yes' 'yes' 'tools/image_processing/imagej2/test-data/analyze_particles_bareoutlines.gif' 'gif' ''
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/analyze_particles_rois.tabular' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Overlay_Masks' 'yes' 'no' 'no' 'tools/image_processing/imagej2/test-data/analyze_particles_overlaymasks.gif' 'gif' 'tools/image_processing/imagej2/test-data/analyze_particles_nothing.tabular'
# Test 5
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/confocal_analyze_particles_rois.tabular' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_threshold_default.tiff' 'yes' '0-Infinity' 0.0 1.0 'Overlay_Masks' 'yes' 'no' 'no' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_threshold_default_overlay_mask.tiff' 'tiff' 'tools/image_processing/imagej2/test-data/confocal_analyze_particles.tabular'

# Skeleton
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'tools/image_processing/imagej2/test-data/skeletonized_blobs.gif' 'no' 'none' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/basic.tabular'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'tools/image_processing/imagej2/test-data/skeletonized_clown.jpg' 'no' 'shortest_branch' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/shortest_branch_basic.tabular'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'tools/image_processing/imagej2/test-data/skeletonized_blobs.gif' 'no' 'none' 'no' 'yes' 'no' 'tools/image_processing/imagej2/test-data/largest_shortest_path_basic.tabular'
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'tools/image_processing/imagej2/test-data/skeletonized_clown.jpg' 'no' 'shortest_branch' 'no' 'yes' 'no' 'tools/image_processing/imagej2/test-data/shortest_branch_all_yes.tabular'

# binary_to_edm
# Test 1
ImageJ --ij2 --headless --debug --jython '/home/ldelisle/Documents/mygit/galaxytools/tools/image_processing/imagej2/imagej2_binary_to_edm_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 1 1 'no' 'no' 'tools/image_processing/imagej2/test-data/blobs_edm.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython '/home/ldelisle/Documents/mygit/galaxytools/tools/image_processing/imagej2/imagej2_binary_to_edm_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 10 3 'yes' 'yes' 'tools/image_processing/imagej2/test-data/blobs_black_edm.gif' 'gif'

# bunwarpj_adapt_transform
# Test 1
bunwarpj -adapt_transform 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/source_elastic_transformation.txt' 'tools/image_processing/imagej2/test-data/adapted_transformation.txt' 2.0

# I don't know why but my new outputs do not pass tests on CI...
# Old do not pass CI neither
# bunwarpj_align
# Test 1
bunwarpj -align 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'NULL' 'tools/image_processing/imagej2/test-data/blobs.gif' 'NULL' 0 2 1 0.0 0.0 1.0 10.0 'tools/image_processing/imagej2/test-data/registered_source1.png' 'tools/image_processing/imagej2/test-data/registered_target1.png'
# Test 2
bunwarpj -align 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'NULL' 'tools/image_processing/imagej2/test-data/blobs.gif' 'NULL' 0 2 1 0.0 0.0 1.0 10.0 'tools/image_processing/imagej2/test-data/registered_source1.png' 'tools/image_processing/imagej2/test-data/registered_target1.png' '-save_transformation'
mv 'tools/image_processing/imagej2/test-data/registered_source1_transf.txt' 'tools/image_processing/imagej2/test-data/source_elastic_transformation_out_full.txt'
mv 'tools/image_processing/imagej2/test-data/registered_target1_transf.txt' 'tools/image_processing/imagej2/test-data/target_elastic_transformation_out_full.txt' 
# Test 3
bunwarpj -align 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/mask_white.png' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/mask_ramp.gif' 0 2 1 0.0 0.0 1.0 10.0 'tools/image_processing/imagej2/test-data/registered_source2.png' 'tools/image_processing/imagej2/test-data/registered_target2.png'

# bunwarpj_compare_elastic_raw
# Test 1
bunwarpj -compare_elastic_raw 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/target_elastic_transformation.txt' 'tools/image_processing/imagej2/test-data/source_raw_transformation.txt' 'tools/image_processing/imagej2/test-data/warping_index_raw_full.txt'  2>&1 | tee 'output_log.txt' && grep -Po 'Warping index = \K[^ ]+' 'output_log.txt' > 'tools/image_processing/imagej2/test-data/warping_index_raw_full.txt'

# bunwarpj_compare_elastic
bunwarpj -compare_elastic 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/target_elastic_transformation.txt' 'tools/image_processing/imagej2/test-data/source_elastic_transformation.txt' 2>&1 | tee 'output_log.txt' && grep -Po 'Warping index = \K[^ ]+' 'output_log.txt' > 'tools/image_processing/imagej2/test-data/warping_index.txt'

# bunwarpj_compare_raw
# Test 1
bunwarpj -compare_raw 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/target_raw_transformation.txt' 'tools/image_processing/imagej2/test-data/source_raw_transformation.txt' 2>&1 | tee 'output_log.txt' && grep -Po 'Warping index = \K[^ ]+' 'output_log.txt' > 'tools/image_processing/imagej2/test-data/warping_index1_full.txt'
# Test 2
bunwarpj -compare_raw 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/source_raw_transformation.txt' 'tools/image_processing/imagej2/test-data/source_raw_transformation.txt' 2>&1 | tee 'output_log.txt' && grep -Po 'Warping index = \K[^ ]+' 'output_log.txt' > 'tools/image_processing/imagej2/test-data/warping_index2.txt'

# bunwarpj_compose_elastic
bunwarpj -compose_elastic 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/target_elastic_transformation.txt' 'tools/image_processing/imagej2/test-data/source_elastic_transformation.txt' 'tools/image_processing/imagej2/test-data/raw_transformation_full.txt'

# bunwarpj_compose_raw_elastic
bunwarpj -compose_raw_elastic 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/target_raw_transformation.txt' 'tools/image_processing/imagej2/test-data/source_elastic_transformation.txt' 'tools/image_processing/imagej2/test-data/composed_raw_elastic_transformation_full.txt'

# bunwarpj_compose_raw
bunwarpj -compose_raw 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/target_raw_transformation.txt' 'tools/image_processing/imagej2/test-data/source_raw_transformation.txt' 'tools/image_processing/imagej2/test-data/composed_raw_transformation_full.txt'

# bunwarpj_convert_to_raw
bunwarpj -convert_to_raw 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/source_elastic_transformation.txt' 'tools/image_processing/imagej2/test-data/source_raw_transformation.txt'

# bunwarpj_elastic_transform
bunwarpj -elastic_transform 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/blobs_direct_transf.txt' 'tools/image_processing/imagej2/test-data/elastic_trans_registered_source1.png'

# bunwarpj_raw_transform
bunwarpj -raw_transform 'tools/image_processing/imagej2/test-data/dotblot.jpg' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/source_raw_transformation.txt' 'tools/image_processing/imagej2/test-data/raw_trans_registered_source1.png'

# Create image
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_create_image_jython_script.py' 'MyTitle' 256 256 1 '8-bit_ramp' 'tools/image_processing/imagej2/test-data/create_image1.jpg'

# Crop
# Test 1
# ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_crop_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 0 0 0 0 1 0 1 0 1 0 'tools/image_processing/imagej2/test-data/blobs.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_crop_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 0 50 0 0 1 0 1 0 1 0 'tools/image_processing/imagej2/test-data/blobs_crop_width50.gif' 'gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_crop_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 0 0 50 0 1 0 1 0 1 0 'tools/image_processing/imagej2/test-data/blobs_crop_top50.gif' 'gif'
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_crop_jython_script.py' 'tools/image_processing/imagej2/test-data/confocal-series-both-channels.tiff' 17 16 18 8 1 1 1 0 1 0 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_cropped.tiff' 'tiff'
# Test 5
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_crop_jython_script.py' 'tools/image_processing/imagej2/test-data/confocal-series-both-channels.tiff' 17 16 18 8 1 0 1 1 1 0 'tools/image_processing/imagej2/test-data/confocal-series-both-channels_cropped_singleZ.tiff' 'tiff'
# Test 6
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_crop_jython_script.py' 'tools/image_processing/imagej2/test-data/confocal-series-both-channels.tiff' 17 16 18 8 1 1 1 1 1 0 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_cropped_singleZ.tiff' 'tiff'

# Enhance contrast
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_enhance_contrast_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'yes' 'None' 'no' 'tools/image_processing/imagej2/test-data/blobs_equalize.gif' 'gif' 
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_enhance_contrast_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' 6.2 'no' 'tools/image_processing/imagej2/test-data/blobs_saturate.gif' 'gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_enhance_contrast_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' 13.0 'yes' 'tools/image_processing/imagej2/test-data/blobs_normalize.gif' 'gif'

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

# Find edges
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_find_edges_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/blobs_find_edges.gif' 'gif'

# Find maxima
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_find_maxima_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'yes' 'no' 10 'Single_Points' 'no' 'no' 'tools/image_processing/imagej2/test-data/blobs_single_points.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_find_maxima_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'yes' 'no' 13 'Maxima_Within_Tolerance' 'no' 'no' 'tools/image_processing/imagej2/test-data/blobs_tolerance.gif' 'gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_find_maxima_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'yes' 'no' 16 'Segmented_Particles' 'yes' 'no' 'tools/image_processing/imagej2/test-data/blobs_segmented.gif' 'gif'
# Test 4
# ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_find_maxima_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'yes' 'no' 10 'List' 'no' 'no' 'tools/image_processing/imagej2/test-data/blobs_list.tabular' 'tabular'
# Test 5
# ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_find_maxima_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'yes' 'no' 10 'Count' 'no' 'no' 'tools/image_processing/imagej2/test-data/blobs_count.tabular' 'tabular'

# make binary
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_make_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/clown.jpg' 1 1 'no' 'no' 'tools/image_processing/imagej2/test-data/clown_binary.jpg' 'jpg'

# Math
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_math_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'Multiply' 'None' 'None' 1.25 'tools/image_processing/imagej2/test-data/blobs_multiply.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_math_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'Min' 'None' 'None' 255.0 'tools/image_processing/imagej2/test-data/blobs_min.gif' 'gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_math_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'Log' 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/blobs_log.gif' 'gif'
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_math_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'Square' 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/blobs_square.gif' 'gif'
# Test 5
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_math_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'Macro' 'v=v+50*sin(d/17)' 'None' 'None' 'tools/image_processing/imagej2/test-data/blobs_macro.gif' 'gif'

# Noise
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_noise_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'gif' 'add_specified_noise' 5.0 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/add_specified_noise.gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_noise_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'gif' 'despeckle' 'None' 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/despeckle.gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_noise_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'gif' 'remove_outliers' 'None' '10.0' '1' 'Bright' 'tools/image_processing/imagej2/test-data/remove_outliers.gif'
# Test 4 fails
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_noise_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif'  'gif' 'remove_nans' 'None' 'None' 'None' 'None' 'tools/image_processing/imagej2/test-data/remove_nans.gif'

# Shadows
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_shadows_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'Northwest' 'tools/image_processing/imagej2/test-data/blobs_northwest.gif' 'gif'

# Sharpen
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_sharpen_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/blobs_sharpen.gif' 'gif'

# Skeletonized 3D
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_skeletonize3d_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' 'tools/image_processing/imagej2/test-data/skeletonized_blobs.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_skeletonize3d_jython_script.py' 'tools/image_processing/imagej2/test-data/clown.jpg' 'no' 'tools/image_processing/imagej2/test-data/skeletonized_clown.jpg' 'jpg'

# Smooth
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_smooth_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'tools/image_processing/imagej2/test-data/blobs_smooth.gif' 'gif'

# Watershed
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_watershed_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' 'tools/image_processing/imagej2/test-data/blobs_watershed_binary.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_watershed_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_threshold_default.tiff' 'yes' 'tools/image_processing/imagej2/test-data/confocal-series-first-channel_default_threshold_watershed.tiff' 'tiff'
