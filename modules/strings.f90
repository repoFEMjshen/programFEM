! parseString
! toUpperCase

      SUBROUTINE parseString(str, del, wordPos, arrSize, wordCount)
      implicit none
      character(len=*), intent(in) :: str
      character, intent(in) :: del
      integer, intent(in) :: arrSize
      integer, intent(out) :: wordPos(arrSize,2)
      integer, intent(out) :: wordCount
      
      integer :: i,wordCountMax
      character :: prevChar, thisChar, del2
      
      ! space and tab alternatively
      del2=del
      if(del.EQ.' ') del2=ACHAR(9)
      if(IACHAR(del).EQ.9) del2=' '
      
      wordCountMax = wordCount
      wordCount = 0
      wordPos(:,:) = 0
      prevChar = del
      do i=1,len(str)
         thisChar = str(i:i)
         if((prevChar.NE.del.AND.prevChar.NE.del2)   &
       .AND.(thisChar.EQ.del.OR.thisChar.EQ.del2)) then
            ! end of a word
            wordPos(wordCount,2) = i - 1
         elseif((prevChar.EQ.del.OR.prevChar.EQ.del2) &
       .AND.(thisChar.NE.del.AND.thisChar.NE.del2))  then
            ! start of a word
            wordCount = wordCount + 1
            wordPos(wordCount,1) = i
         endif
         prevChar = thisChar
      enddo
      if(prevChar.NE.del.AND.prevChar.NE.del2) then
         !end of string ends the word
         wordPos(wordCount,2) = len(str)
      endif
      END SUBROUTINE
	  
      
      SUBROUTINE toUpperCase(ch)
      implicit none
      character, intent(inout) :: ch
      integer :: ascii
      ascii = IACHAR(ch)
      if(ascii.GE.97.AND.ascii.LE.122) ch=ACHAR(ascii-32)
      END SUBROUTINE
