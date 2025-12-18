      complex(8) function zsign(a,b)
      implicit none
      complex(8) a,b
      if(dble(a)*dble(b)>0)then
         zsign=a
      else
         zsign=-a
      endif
      return
      end
