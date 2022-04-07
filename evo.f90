       module evo
         integer parameter maxev=100
         integer evl(maxev,2)
         integer currev
         subroutine putev(target,time,*)
         integer target,time
         integer i
         logical done
         done=false
         do i=1,maxev
            if (evl(i,1)=0) then
               evl(i,1)=target
               evl(i,2)=time
               done=.true.
            end if
         end do
         if (done=.false.) then
            
         end if 
         
         end subroutine putev
       end module evo