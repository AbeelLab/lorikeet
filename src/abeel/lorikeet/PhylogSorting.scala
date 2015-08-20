package abeel.lorikeet

import atk.compbio.tree.Tree
import atk.util.Tool
import scala.collection.JavaConversions._
import java.io.PrintWriter
import java.io.File

object PhylogSorting extends Tool {

  override val version = """
            2014/12/16  Initial version
                        Improved handling of missing lineages
                        Exposed agreement fraction and SNP distance parameters for users
    
    
  """

  case class Config(input: File = null, output: File = null, tree: File = null, snpmatrix: File = null, maxDist: Int = 500, consensusFraction: Double = 0.6)

  def main(args: Array[String]): Unit = {
    val parser = new scopt.OptionParser[Config]("java -jar lorikeet.jar fix-lineages") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input lineage information. (Output of merge-spoligotypes)")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output file.")
      opt[File]('t', "tree") required () action { (x, c) => c.copy(tree = x) } text ("Phylogenetic tree file in NWK format.")
      opt[File]('s', "snpmatrix") required () action { (x, c) => c.copy(snpmatrix = x) } text ("Matrix with pairwise SNP distances (created with pelican)")
      opt[Int]("distance") action { (x, c) => c.copy(maxDist = x) } text ("Maximum distance to consider closest neighbors. [Default=500]")
      opt[Double]("fraction") action { (x, c) => c.copy(consensusFraction = x) } text ("Fraction of closest neighbors that need to agree to perform change. [Default=0.6]")
    }

    parser.parse(args, Config()).map { config =>
      fix(config);
    }

  }

  private def fix(config: Config) {
    val tree = new Tree(config.tree.toString())

    val output = "correctedlineages.txt"
    val lines = tMap(tLines(config.input), 0, 3)
    val lines2 = tMap(tLines(config.input))

    val rowMap = tMap(tLines(config.snpmatrix))
    val header = rowMap("Taxon").split("\t")

    val countPerLineage = lines.groupBy(_._2).mapValues(_.size)
    val pw = new PrintWriter(config.output)
    pw.println(generatorInfo)
    val leaves = tree.getLeaves(tree.root)
    var count = 0
    var unsupported = 0
    for (l <- leaves) {
      val name = l.getName()
      // find 10 closest
      val list = header.zip(rowMap(l.getName).split("\t").map(_.toInt))
      val closest10 = list.sortBy(_._2).drop(1).filter(f => lines.getOrElse(f._1, "MISSING").startsWith("LIN")).filter(f => f._2 < config.maxDist).take(10)

      // average distance of closests

      // proposed renaming, if any

      val x1 = closest10.map(f => lines(f._1) -> f._2)
      val x2 = x1.groupBy(_._1)
      //      println(x2)
      val x3 = x2.mapValues(f => f.size -> (f.map(_._2).sum / f.size.toDouble))

      val consensus = x3.filter(_._2._1 > config.consensusFraction * closest10.size)

      /**
       * Lineage not supported by neighbors
       *
       */
      val lin = lines.getOrElse(name, "MISSING\tMISSING\tMISSING")
      val lin2 = lines2.getOrElse(name, "MISSING\tMISSING\tMISSING")
      if (x3.size >= 1 && consensus.size == 0 && countPerLineage(lin) >= closest10.size) {

        println("Unsupported: " + name + "\t" + lin + "\t" + x3 + "\t" + countPerLineage(lin))
        pw.println(name + "\t" + lines2(name).split("\t").take(2).mkString("\t") + "\t-")
        unsupported += 1
      } else if (consensus.size >= 1 && consensus.head._2._2 < config.maxDist && !consensus.head._1.equals(lin)) {

        count += 1
        println("Rename: " + name + "\t" + lin + "\t" + x3 + "\t->" + consensus.head._1)
        pw.println(name + "\t" + lin2.split("\t").take(2).mkString("\t") + "\t" + consensus.head._1)
      } else {
        pw.println(name + "\t" + lin2)
      }

    }
    println("Corrections: " + count)
    println("Unsupported: " + unsupported)

    pw.close

  }

}