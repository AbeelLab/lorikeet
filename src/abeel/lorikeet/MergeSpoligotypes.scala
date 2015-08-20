package abeel.lorikeet

import java.io.File
import java.io.PrintWriter
import atk.io.PatternFileFilter
import atk.util.Tool
import atk.compbio.tree.Tree
import atk.io.IOTools

object MergeSpoligotypes extends Tool {

  override val version = """
    Program to merge many spoligotype results and and infer lineage information
    
    2014/11/21  Initial version with version information
                Added recursive option
                Use 43 markers instead of 104

    2015/08/20  Initial public version
    
    
   """

  def main(args: Array[String]): Unit = {
    case class Config(val input: Seq[File] = Seq(),
      val output: String = null,
      val pattern: String = """.*.spoligotype""",
      val recursive: Boolean = false)
    val spacersToUse = 43

    val parser = new scopt.OptionParser[Config]("java -jar lorikeet.jar merge-spoligotypes") {
      opt[File]('i', "input") required () unbounded () action { (x, c) => c.copy(input = c.input :+ x) } text ("Input directory that contains all spoligotype files. You can specify multiple -i arguments")

      opt[String]('o', "output") action { (x, c) => c.copy(output = x) } text ("Output prefix")
      opt[Unit]('r', "recursive") action { (x, c) => c.copy(recursive = true) } text ("Search input directories recursively [Default=true]")
      opt[String]('p', "pattern") action { (x, c) => c.copy(pattern = x) } text ("File name pattern for the input files. [Default=\".*.spoligotype]\"")

    }
    parser.parse(args, Config()) map { config =>
      val spoligoList = config.input.map { f =>

        if (config.recursive)
          IOTools.recurse(f, new PatternFileFilter(config.pattern)).toList
        else
          f.listFiles(new PatternFileFilter(config.pattern)).toList
      }.flatten.toList
      println(spoligoList)
      condense(config.output, spoligoList)

      val pw = new PrintWriter(config.output + "spoligotype.matrix")
      pw.println(generatorInfo)
      val str = new Range(1, spacersToUse + 1, 1).toList.mkString("\t")
      pw.println("$$\t" + str)

      val bytes = "1111111111 1111111110 0111111111 1100001111 111".replaceAll(" ", "").getBytes().take(spacersToUse)
      //      bytes(strains.indexOf(p._2)) = '1'
      pw.println("H37RV\t" + bytes.map(_.toChar).mkString("\t"))

      spoligoList.map(p => {
        val content = tLines(p).takeRight(3)(0).replaceAll(" ", "")

        val bytes = content.getBytes().take(spacersToUse)
        pw.println(p.getName + "\t" + bytes.map(_.toChar).mkString("\t"))

      })
      pw.close

    }

  }
  def condense(outputPrefix: String, listOfFiles: List[File]) = {

    val pw = new PrintWriter(outputPrefix + ".spoligotypes.txt");
    val pw2 = new PrintWriter(outputPrefix + ".lineages.txt");
    pw.println(generatorInfo)
    pw2.println(generatorInfo)

    for (file <- listOfFiles) {
      val lines = tLines(file)
      println("Processing: " + file)
      pw.println(file.getName() + "\t" + lines.takeRight(2).take(1)(0))
      val spoligo = lines.takeRight(1)(0)
      val lin = Spoligotype2Lineage.convert(spoligo.split("\t")(0).split(",")(0))

      pw2.println(file.getName() + "\t" + spoligo + "\t" + lin)

    }
    pw.close
    pw2.close
  }
}