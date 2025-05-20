# Rust implementation of the TCS pipeline

Rust implementation with parallel processing, performance improvments and more functionalilties

## Usage

```
Usage: tcs <COMMAND>

Commands:
  run        Run the TCS pipeline
  generate   Generate a param file through CLI
  dr         Run the TCS HIV-1 DR Pipeline,
  dr-params  List param for the DR pipeline, w/o aurguments it will list all available version numbers
  sdrm       SDRM pipeline followed by HIV-1 DR pipeline
  log        Aggregate log files and reorganize the directory structure after TCS or DR pipeline
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

```

### Run the TCS pipeline

```
Usage: tcs run [OPTIONS] --input <INPUT> --param <PARAM>

Options:
-i, --input <INPUT> Input directory path
-p, --param <PARAM> param file path
--keep-original keep original files
-h, --help Print help
```

### Run the TCS_DR pipeline

```
Usage: tcs dr [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>      Input directory path
  -v, --version <VERSION>  DR version number [default: v1]
      --keep-original      keep original files
  -h, --help               Print help
```

### SDRM pipeline followed by HIV-1 DR pipeline

```
Usage: tcs sdrm --input <INPUT>

Options:
  -i, --input <INPUT>      Input directory path
  -v, --version <VERSION>  DR version number [default: v1]
  -h, --help               Print help
```

### Aggregate log files and reorganize the directory structure after TCS or DR pipeline

```
Usage: tcs log --input <INPUT>

Options:
  -i, --input <INPUT>  Input directory path
  -h, --help           Print help
```

### List param for the DR pipeline, w/o aurguments it will list all available version numbers

```
Usage: tcs dr-params

Options:
  -v, --version <VERSION>  Print out params for a specific version
  -h, --help               Print help
```

### Generate a param file through CLI

```
Usage: tcs generate

Options:
  -h, --help  Print help
```
