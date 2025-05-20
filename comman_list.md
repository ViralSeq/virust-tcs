# TCS parameter generator CLI

`tcs generator` or `tcs g`

# TCS DR pipeline

`tcs dr -v <dr_version> -i <input_dir>`
`tcs dr -v <dr_version> -i <input_dir> --keep-original`

## addtional commands

`tcs dr list` #show all availabe DR pipelines.
`tcs dr list -v <dr_version>` # print out params for the specific version

# TCS general pipeline

`tcs run -p <param.json> -i <input_dir>`
`tcs run -p <param.json> -i <input_dir> --keep-original`

# TCS log pipeline

`tcs log -i <dir_to_batch_tcs_jobs>`

# TCS hiv sdrm pipeline

`tcs sdrm -v <dr_version> -i <input_dir>`

# general

`tcs -h` # help function
`tcs -v` # print out version
