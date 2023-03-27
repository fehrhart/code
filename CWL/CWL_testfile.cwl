#Example on a data mangling script, data input is defined in the script

cwlVersion: v1.0
class: CommandLineTool

requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: GeneHancer.R 
      entry:
        $include: .\GeneHancer.R

baseCommand: ["Rscript", "GeneHancer.R"]

inputs: []
#input defined in the analysis script, so not needed here
 
outputs:
  output_file:
  type: File
  outputBinding:
   glob: GeneHancerScore5_input.txt
  
#Steps:
#only needed if there is more than one script, otherwise, steps are defined CWL files
