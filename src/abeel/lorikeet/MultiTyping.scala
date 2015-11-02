package abeel.lorikeet

import java.io.File
import java.io.PrintWriter
import atk.io.PatternFileFilter
import atk.util.Tool
import atk.compbio.tree.Tree
import atk.io.IOTools

object MultiTyping extends Tool {

  override val version = """
    Program to merge many spoligotype results, condense multiple libraries and infer lineage information
    
    This is a conceptual successor to 'merge-spoligotypes'
    
    2015/10/30  Initial version with version information
                
    
    
   """
  case class Config(val input: Seq[File] = Seq(),
    val output: String = null,
    val pattern: String = """.*.spoligotype""",
    val recursive: Boolean = false)
  def main(args: Array[String]): Unit = {

    val spacersToUse = 43

    val parser = new scopt.OptionParser[Config]("java -jar lorikeet.jar multi-typing") {
      opt[File]('i', "input") required () unbounded () action { (x, c) => c.copy(input = c.input :+ x) } text ("Input directory that contains all spoligotype files. You can specify multiple -i arguments")

      opt[String]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output prefix")
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
      condense(config.output, spoligoList, config)

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
  def condense(outputPrefix: String, listOfFiles: List[File], config: Config) = {

    val pw = new PrintWriter(outputPrefix + ".spoligotypes.txt");
//    val pw2 = new PrintWriter(outputPrefix + ".lineages.txt");
    pw.println(generatorInfo)
//    pw2.println(generatorInfo)

    for (file <- listOfFiles) {
      val lines = tLines(file)
      println("Processing: " + file)
      println(file.getName() + "\t" + lines.takeRight(2).take(1)(0))
      //      val spoligo = lines.takeRight(1)(0)
      val groups = lines.grouped(46).toList
      println(groups.size)

      println(groups.map(_(45)))

      val summed = (groups.map(_.take(43).map(l => l.split("\t")(1).toInt))).transpose.map(_.sum)
      val sorted = summed.sorted
      // threshold is 5% of 3 highest values
      val threshold = math.max(10, (sorted(40) + sorted(41) + sorted(42)) / 60)
      println(threshold)
      val binary = summed.map(f => if (f >= threshold) 1 else 0)

      //      val lin = Spoligotype2Lineage.convert(spoligo.split("\t")(0).split(",")(0))
      println(summed.mkString(", "))
      //      pw2.println(file.getName() + "\t" + spoligo + "\t" + lin)

      val buffer = binary.grouped(10).toList.map(f => f.mkString("")).mkString(" ")
      //     val binStr=buffer
      //val buffer=new StringBuffer

     println("## Binary spoligotype: \n" + buffer.toString())
      //      val binary = buffer.toString.replaceAll(" ", "").toList.take(43).map(f => f.toString.toInt)
      val octal = binary.grouped(3).map(oc => if (oc.length == 3) (oc(0) * 4 + oc(1) * 2 + oc(2) * 1) else oc(0)).toList
      println("## Octal spoligotype: \n" + octal.mkString(""))
      println("## Lineage:\n" + Lineages.lineage(buffer.toString))
      println()

//      pw.println(new File(config.input))

      val parent=file.getParentFile().getCanonicalFile().getName()
      val grantparent=file.getParentFile().getCanonicalFile().getParentFile().getCanonicalFile().getName()
      
      pw.println(grantparent+"\t"+parent+"\t"+buffer.toString()+"\t"+octal.mkString("")+"\t"+Lineages.lineage(buffer.toString)+"\t"+summed.mkString(","))
      
    }
    pw.close
//    pw2.close
  }
}