/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef.FilePath
import htsjdk.samtools.util.{CoordMath, OverlapDetector}

/** Convenience methods for [[Amplicon]] */
object Amplicon {
  /** Builds a [[OverlapDetector]] for the given amplicons. */
  def overlapDetector(amplicons: Iterator[Amplicon]): OverlapDetector[Amplicon] = {
    val detector = new OverlapDetector[Amplicon](0,0)
    amplicons.foreach(amp => detector.addLhs(amp, amp))
    detector
  }

  /** Builds a [[OverlapDetector]] for the given file of amplicons. */
  def overlapDetector(path: FilePath): OverlapDetector[Amplicon] = {
    Amplicon.overlapDetector(amplicons=Metric.iterator[Amplicon](path))
  }
}


// Developer note: the left/right start/end members are snakecase, which is the convention for  "Metric" files
// (originally required by [[TrimPrimers]]).  The are made private so that code that uses [[Amplicon]] uses the
// camel case instead (leftStart/leftEnd/rightStart/rightEnd).
/** A Locatable Amplicon class.
  * @param chrom the chromosome for the amplicon
  * @param left_start the 1-based start position of the left-most primer
  * @param left_end the 1-based end position inclusive of the left-most primer
  * @param right_start the 1-based start position of the right-most primer
  * @param right_end the 1-based end position inclusive of the right-most primer
  */
case class Amplicon
( chrom: String,
  private val left_start: Option[Int],
  private val left_end: Option[Int],
  private val right_start: Option[Int],
  private val right_end: Option[Int],
  private val id: Option[String] = None
) extends GenomicSpan with Metric {
  @inline def contig: String  = chrom

  val (s, e, lpl) = (left_start, left_end, right_start, right_end) match {
    case (Some(left_start), Some(left_end), Some(right_start), Some(right_end)) =>
      require(left_start <= left_end, f"leftStart is > leftEnd: $this")
      require(right_start <= right_end, f"rightStart is > rightEnd: $this")
      require(left_start <= right_start, f"leftStart is > rightStart: $this")

      def leftPrimerLength: Int       = CoordMath.getLength(left_start, left_end)
      def rightPrimerLength: Int      = CoordMath.getLength(right_start, right_end)
      def longestPrimerLength: Int    = Math.max(leftPrimerLength, rightPrimerLength)
      def leftPrimerLocation: String  = f"$chrom:$left_start-$left_end"
      def rightPrimerLocation: String = f"$chrom:$right_start-$right_end"
      def ampliconLocation: String    = f"$chrom:$left_start-$right_end"
      def identifier: String          = this.id.getOrElse(ampliconLocation)
      (left_start, right_end, Math.max(leftPrimerLength, rightPrimerLength))
    case (Some(left_start), Some(left_end), None, None) =>
      require(left_start <= left_end, f"leftStart is > leftEnd: $this")
      @inline def leftStart: Int  = left_start
      @inline def leftEnd: Int    = left_end
      def leftPrimerLength: Int       = CoordMath.getLength(leftStart, leftEnd)
      def rightPrimerLength: Int      = 0
      def longestPrimerLength: Int    = leftPrimerLength
      def leftPrimerLocation: String  = f"$chrom:$leftStart-$leftEnd"
      (left_start, leftEnd, leftPrimerLength)
    case (None, None, Some(right_start), Some(right_end)) => 
      require(right_start <= right_end, f"rightStart is > rightEnd: $this")
      @inline def rightStart: Int = right_start
      @inline def rightEnd: Int   = right_end

      def leftPrimerLength: Int       = 0
      def rightPrimerLength: Int      = CoordMath.getLength(rightStart, rightEnd)
      def rightPrimerLocation: String = f"$chrom:$rightStart-$rightEnd"
      (right_start, rightEnd, rightPrimerLength)
    case _ =>
      throw new Exception(s"At least (left_start and left_end) or (right_start and right_end) need to be set in every row of the primer file.")
  }
  @inline def longestPrimerLength: Int  = lpl
  @inline def start: Int = s
  @inline def end: Int = e
}
