remotefiles_dir='flatMISMIP_testbeds/long_models_yang/'
localfiles_dir='../long_models_yang/'
folder1='model_W5000_GL400_FC30000/'
folder2='model_W5000_GL400_FC120000/'
folder3='model_W11000_GL400_FC30000/'
folder4='model_W11000_GL400_FC120000/'
file_wildcard='*Local*Extended.mat'
sftp mycluster <<EOF
# get $remotefiles_dir$folder1$file_wildcard $localfiles_dir$folder1
get $remotefiles_dir$folder2$file_wildcard $localfiles_dir$folder2
get $remotefiles_dir$folder3$file_wildcard $localfiles_dir$folder3
get $remotefiles_dir$folder4$file_wildcard $localfiles_dir$folder4
EOF  
