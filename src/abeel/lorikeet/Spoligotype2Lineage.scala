package abeel.lorikeet

import atk.util.Tool
import scala.io.Source

object Spoligotype2Lineage extends Tool {
  private val stream = Spoligotype2Lineage.getClass().getResourceAsStream("/spoligotype2lineage.txt")
  private val lines = Source.fromInputStream(stream).getLines

  private val map = tMap(lines.filterNot(f => f.startsWith("#")).filterNot(f => f.trim.size == 0).toList)
  println(lines.mkString("\n"))
  def convert(str: String): String = {
    map.getOrElse(str, str)

  }

}

