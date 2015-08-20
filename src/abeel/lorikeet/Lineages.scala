package abeel.lorikeet

import scala.io.Source
import java.io.File
import java.io.PrintWriter
import scala.collection.immutable.BitSet
import scala.collection.mutable.HashSet
import atk.util.Tool

/*
 * 
 * 
 */
object Lineages extends Tool {

  class Spoligotype(key: String) {

    lazy val string = key.replaceAll(" ", "")

    lazy val binary = BitSet(string.zipWithIndex.filter(p => p._1 == '1').map(_._2): _*)

  }

  /**
   * mask is the absence and presence of markers
   *
   * used is a bitset indicating which markers need to be used
   */
  class SpoligoTypeLineage(val name: String, val mask: BitSet, val sample: String) {
    val used = BitSet((0 to 42): _*)
    def distance(spol: Spoligotype): Int = {
      val maskResult = (spol.binary ^ mask)
      val result = maskResult & used
      result.size
    }
  }
  val input = Source.fromInputStream(Lineages.getClass().getResourceAsStream("/mytypes.txt"))("iso-8859-1").getLines().toList.filterNot(l=>l.startsWith("#")||l.trim.size==0)
  val lineages = input.map(f => {
    val arr = f.split("[\\s]+")
    val name = arr(0)
    new SpoligoTypeLineage(name, BitSet(arr(1).zipWithIndex.filter(p => p._1 == 'n').map(_._2): _*), f)
  })


  def lineage(inputLine: String): String = {

  
    val spol = new Spoligotype(inputLine);
    val distances = lineages.map(f => (f, f.distance(spol))).sortBy(_._2)

    val smallestDistance = distances.head._2

    val ll = distances.filterNot(f => f._2 > smallestDistance || f._2>2).map(_._1.name).toSet
   
    (if(ll.size==0)"UNKNOWN" else ll.mkString(",")) + "\t" + distances.toList.sortBy(_._2).map(f => f._1.name + " (" + f._2 + ")").take(10).toSet.mkString(";")
   

  }

}