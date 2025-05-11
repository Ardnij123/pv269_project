version 1.0



task plot_compare {
  input {
    Array[File] fasta_files
    Int resolution = 100
    String mdp_options = ""
    Boolean only_compare = true
  }

  String comp = if only_compare then "--compare-only" else "--compare"

  command <<<
  moddotplot static -f ~{sep=' ' fasta_files} --compare-only ~{mdp_options} -r ~{resolution} -o "similiarity-plots" ~{comp} --width 1 --dpi ~{resolution * 5}
  >>>

  runtime {
    docker: "jindmen/moddotplot:latest"
    disks: "local-disk 20 SSD"
    cpu: 2
    memory: "2GB"
    preemptible: 1
  }

  output {
    Array[File] sim_files = glob("similiarity-plots/*")
    Array[File] comp_files = glob("similiarity-plots/*_COMPARE.png")
  }
}



task split_fasta {
  input {
    File fasta_file
  }

  command <<<
    seqkit split --by-id --out-dir 'fasta-splitted' --threads 4 '~{fasta_file}'
  >>>

  runtime {
    docker: "staphb/seqkit:latest"
    disks: "local-disk 20 SSD"
    cpu: 6
    memory: "8 GB"
    preemptible: 1
  }

  output {
    Array[File] split_files = glob("fasta-splitted/*")
  }
}



task grid_output {
  input {
    File fasta_file
    Array[File] plots
    Int resolution = 500
  }

  command <<<
    sequences=$(\
      grep ">" ~{fasta_file} |
      sed 's/^.//' |
      sed 's/chr//' | sed 's/^X/9999X/' | sed 's/^Y/9999Y/'|
      sort -n |
      sed 's/^9999//' | sed 's/^/chr/'
    )

    mkdir "temp-plots"

    echo "$sequences" | while read seq; do
      echo "$sequences" | while read seq2; do
        file="`grep "${seq}_${seq2}_COMPARE.png" ~{write_lines(plots)}`"
        target="temp-plots/${seq}_${seq2}.png"
        if [ ! "$file" ]; then
          convert -size "~{resolution}x~{resolution}" xc:transparent "$target"
        else
          ln -s "$file" "$target"
        fi
      done
    done

    echo "$sequences" | while read seq; do
      convert -append `echo "$sequences" | while read seq2; echo "temp-plots/${seq}_${seq2}.png"; done` "row_$seq.png"
    done

    convert +append `echo "$sequences" | while read seq; do echo "row_$seq.png"; done | tac` "grid.png"
  >>>

  runtime {
    docker: "dpokidov/imagemagick:latest-bookworm"
    disks: "local-disk 20 SSD"
    cpu: 4
    memory: "8 GB"
    preemptible: 1
  }

  output {
    File grid_file = "grid.png"
  }
}



workflow moddotplot_cross {
  input {
    File fasta_file
    Boolean create_grid = true
    String mdp_options = ""
    Int resolution = 100
  }

  call split_fasta {
    input:
    fasta_file = fasta_file
  }

  scatter ( seqs in cross(split_fasta.split_files, split_fasta.split_files) )
  {
    call plot_compare {
      input:
      fasta_files = [ seqs.left, seqs.right ],
      resolution = resolution,
      mdp_options = "~{mdp_options} --no-bedpe",
      only_compare = true
    }
  }

  if (create_grid) {
    call grid_output {
      input:
      fasta_file = fasta_file,
      plots = flatten(plot_compare.comp_files),
      resolution = resolution * 5
    }
  }

  output {
    Array[File] sim_plots = flatten(plot_compare.sim_files)
    File? grid_file = grid_output.grid_file
  }
}
