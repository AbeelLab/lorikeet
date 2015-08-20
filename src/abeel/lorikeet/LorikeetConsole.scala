package abeel.lorikeet

object LorikeetConsole {

  def main(args: Array[String]): Unit = {

    if (args.length == 0) {

      listInstructions
    } else {
      args(0) match {
        case "list" => listInstructions
        case "help" => listInstructions

        case "spoligotype" => LorikeetSpoligotyper.main(args.drop(1))
        case "merge-spoligotypes" => MergeSpoligotypes.main(args.drop(1))
        case "fix-lineages" => PhylogSorting.main(args.drop(1))
        case "_" => listInstructions
      }
    }

  }

  def listInstructions() {
    println("Usage:java -jar lorikeet.jar [instruction] [instruction options...]")
    println("Instructions:")
    println("\tspoligotype            Spoligotype BAM file based on digital spoligotyping")
    println("\tmerge-spoligotypes     Merge multiple spoligotype files together in a single file")
    println("\tfix-lineages           Fix lineage annotations based on phylogenetic tree and SNP distance matrix using a KNN classifier. ")

  }

}