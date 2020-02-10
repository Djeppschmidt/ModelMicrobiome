
#define species functions
# 100 high response, linear functions Factor 1 ####
# 50 positive, 50 negative

#linear high responses factor A
spp1<-function(a,b,c,d,e) {(c*d*e)-(a-b)^2 +(0*(a+b+c+d+e))}
spp2<-function(a,b,c,d,e) {(c*d*e)-(a-c)^2+(0*(a+b+c+d+e))}
spp3<-function(a,b,c,d,e) {(c*d*e)-(a-d)^2+(0*(a+b+c+d+e))}
spp4<-function(a,b,c,d,e) {(c*d*e)-(a-e)^2+(0*(a+b+c+d+e))}
spp5<-function(a,b,c,d,e) {(c*d*e)-(a+b)^2+(0*(a+b+c+d+e))}
spp6<-function(a,b,c,d,e) {(c*d*e)-(a+c)^2+(0*(a+b+c+d+e))}
spp7<-function(a,b,c,d,e) {(c*d*e)-(a+d)^2+(0*(a+b+c+d+e))}
spp8<-function(a,b,c,d,e) {(c*d*e)-(a+e)^2+(0*(a+b+c+d+e))}
spp9<-function(a,b,c,d,e) {(c*d*e)-(b)^2+(0*(a+b+c+d+e))}
spp10<-function(a,b,c,d,e) {(c*d*e)-(a)^2+(0*(a+b+c+d+e))}

spp11<-function(a,b,c,d,e) {(d*e)-(a-b)^2 +(0*(a+b+c+d+e))}
spp12<-function(a,b,c,d,e) {(d*e)-(a-c)^2+(0*(a+b+c+d+e))}
spp13<-function(a,b,c,d,e) {(d*e)-(a-d)^2+(0*(a+b+c+d+e))}
spp14<-function(a,b,c,d,e) {(d*e)-(a-e)^2+(0*(a+b+c+d+e))}
spp15<-function(a,b,c,d,e) {(d*e)-(a+b)^2+(0*(a+b+c+d+e))}
spp16<-function(a,b,c,d,e) {(d*e)-(a+c)^2+(0*(a+b+c+d+e))}
spp17<-function(a,b,c,d,e) {(d*e)-(a+d)^2+(0*(a+b+c+d+e))}
spp18<-function(a,b,c,d,e) {(b*d*e)-(a+e)^2+(0*(a+b+c+d+e))}
spp19<-function(a,b,c,d,e) {(d*e)-(b)^2+(0*(a+b+c+d+e))}
spp20<-function(a,b,c,d,e) {(d*e)-(a)^2+(0*(a+b+c+d+e))}

spp21<-function(a,b,c,d,e) {(a*b)-(a-b)^2 +(0*(a+b+c+d+e))}
spp22<-function(a,b,c,d,e) {(a*b)-(a-c)^2+(0*(a+b+c+d+e))}
spp23<-function(a,b,c,d,e) {(a*b)-(a-d)^2+(0*(a+b+c+d+e))}
spp24<-function(a,b,c,d,e) {(a*b)-(a-e)^2+(0*(a+b+c+d+e))}
spp25<-function(a,b,c,d,e) {(a*b*c)-(a+b)^2+(0*(a+b+c+d+e))}
spp26<-function(a,b,c,d,e) {(a*b)-(a+c)^2+(0*(a+b+c+d+e))}
spp27<-function(a,b,c,d,e) {(a*b)-(a+d)^2+(0*(a+b+c+d+e))}
spp28<-function(a,b,c,d,e) {(a*b)-(a+e)^2+(0*(a+b+c+d+e))}
spp29<-function(a,b,c,d,e) {(a*b)-(b)^2+(0*(a+b+c+d+e))}
spp30<-function(a,b,c,d,e) {(a*b)-(a)^2+(0*(a+b+c+d+e))}

spp31<-function(a,b,c,d,e) {(c*e)-(a-b)^2 +(0*(a+b+c+d+e))}
spp32<-function(a,b,c,d,e) {(c*e)-(a-c)^2+(0*(a+b+c+d+e))}
spp33<-function(a,b,c,d,e) {(c*e)-(a-d)^2+(0*(a+b+c+d+e))}
spp34<-function(a,b,c,d,e) {(c*e)-(a-e)^2+(0*(a+b+c+d+e))}
spp35<-function(a,b,c,d,e) {(c*e)-(a+b)^2+(0*(a+b+c+d+e))}
spp36<-function(a,b,c,d,e) {(c*e)-(a+c)^2+(0*(a+b+c+d+e))}
spp37<-function(a,b,c,d,e) {(c*e)-(a+d)^2+(0*(a+b+c+d+e))}
spp38<-function(a,b,c,d,e) {(c*d*e)-(a+e)^2+(0*(a+b+c+d+e))}
spp39<-function(a,b,c,d,e) {(c*e)-(b)^2+(0*(a+b+c+d+e))}
spp40<-function(a,b,c,d,e) {(c*e)-(a)^2+(0*(a+b+c+d+e))}

spp41<-function(a,b,c,d,e) {(c*d)-(a-b)^2 +(0*(a+b+c+d+e))}
spp42<-function(a,b,c,d,e) {(c*d)-(a-c)^2+(0*(a+b+c+d+e))}
spp43<-function(a,b,c,d,e) {(c*d)-(a-d)^2+(0*(a+b+c+d+e))}
spp44<-function(a,b,c,d,e) {(c*d)-(a-e)^2+(0*(a+b+c+d+e))}
spp45<-function(a,b,c,d,e) {(c*d)-(a/d+b/d)^2+(0*(a+b+c+d+e))}
spp46<-function(a,b,c,d,e) {(c*d*e)-(a+c)^2+(0*(a+b+c+d+e))}
spp47<-function(a,b,c,d,e) {(c*d)-(a+d)^2+(0*(a+b+c+d+e))}
spp48<-function(a,b,c,d,e) {(c*d)-(a+e/b)^2+(0*(a+b+c+d+e))}
spp49<-function(a,b,c,d,e) {(c*d)-(b)^2+(0*(a+b+c+d+e))}
spp50<-function(a,b,c,d,e) {(c*d)-(a)^2+(0*(a+b+c+d+e))}

spp51<-function(a,b,c,d,e) {(b*e)-(a-b)^2 +(0*(a+b+c+d+e))}
spp52<-function(a,b,c,d,e) {(b*e)-(a-c)^2+(0*(a+b+c+d+e))}
spp53<-function(a,b,c,d,e) {(b*e)-(a-d)^2+(0*(a+b+c+d+e))}
spp54<-function(a,b,c,d,e) {(b*e)-(a-e)^2+(0*(a+b+c+d+e))}
spp55<-function(a,b,c,d,e) {(b*e)-(a+b)^2+(0*(a+b+c+d+e))}
spp56<-function(a,b,c,d,e) {(b*e)-(a+c)^2+(0*(a+b+c+d+e))}
spp57<-function(a,b,c,d,e) {(b*e)-(a+d)^2+(0*(a+b+c+d+e))}
spp58<-function(a,b,c,d,e) {(b*e)-(a+e)^2+(0*(a+b+c+d+e))}
spp59<-function(a,b,c,d,e) {(b*e)-(b)^2+(0*(a+b+c+d+e))}
spp60<-function(a,b,c,d,e) {(b*e)-(a)^2+(0*(a+b+c+d+e))}

spp61<-function(a,b,c,d,e) {(a*e)-(a-b)^2 +(0*(a+b+c+d+e))}
spp62<-function(a,b,c,d,e) {(a*e)-(a-c)^2+(0*(a+b+c+d+e))}
spp63<-function(a,b,c,d,e) {(a*e)-(a-d)^2+(0*(a+b+c+d+e))}
spp64<-function(a,b,c,d,e) {(a*e)-(a-e)^2+(0*(a+b+c+d+e))}
spp65<-function(a,b,c,d,e) {(a*d*e)-(a+b)^2+(0*(a+b+c+d+e))}
spp66<-function(a,b,c,d,e) {(a*e)-(a+c)^2+(0*(a+b+c+d+e))}
spp67<-function(a,b,c,d,e) {(a*e)-(a+d)^2+(0*(a+b+c+d+e))}
spp68<-function(a,b,c,d,e) {(a*c*e)-(a+e)^2+(0*(a+b+c+d+e))}
spp69<-function(a,b,c,d,e) {(a*e)-(b)^2+(0*(a+b+c+d+e))}
spp70<-function(a,b,c,d,e) {(a*e)-(a)^2+(0*(a+b+c+d+e))}

spp71<-function(a,b,c,d,e) {(a*c)-(a-b)^2 +(0*(a+b+c+d+e))}
spp72<-function(a,b,c,d,e) {(a*c)-(a-c)^2+(0*(a+b+c+d+e))}
spp73<-function(a,b,c,d,e) {(a*c)-(a-d)^2+(0*(a+b+c+d+e))}
spp74<-function(a,b,c,d,e) {(a*c)-(a-e)^2+(0*(a+b+c+d+e))}
spp75<-function(a,b,c,d,e) {(a*c)-(a+b/c)^2+(0*(a+b+c+d+e))}
spp76<-function(a,b,c,d,e) {(a*c)-(a+c/c)^2+(0*(a+b+c+d+e))}
spp77<-function(a,b,c,d,e) {(a*c)-(a+d/c)^2+(0*(a+b+c+d+e))}
spp78<-function(a,b,c,d,e) {(a*c)-(a+e/c)^2+(0*(a+b+c+d+e))}
spp79<-function(a,b,c,d,e) {(a*c)-(b)^2+(0*(a+b+c+d+e))}
spp80<-function(a,b,c,d,e) {(a*c)-(a)^2+(0*(a+b+c+d+e))}

spp81<-function(a,b,c,d,e) {(b*c)-(a-b)^2 +(0*(a+b+c+d+e))}
spp82<-function(a,b,c,d,e) {(b*c)-(a-c)^2+(0*(a+b+c+d+e))}
spp83<-function(a,b,c,d,e) {(b*c)-(a-d)^2+(0*(a+b+c+d+e))}
spp84<-function(a,b,c,d,e) {(b*c)-(a-e)^2+(0*(a+b+c+d+e))}
spp85<-function(a,b,c,d,e) {(b*c)-(a+b/a)^2+(0*(a+b+c+d+e))}
spp86<-function(a,b,c,d,e) {(b*c)-(a+c)^2+(0*(a+b+c+d+e))}
spp87<-function(a,b,c,d,e) {(b*c)-(a+d)^2+(0*(a+b+c+d+e))}
spp88<-function(a,b,c,d,e) {(b*c)-(a+e)^2+(0*(a+b+c+d+e))}
spp89<-function(a,b,c,d,e) {(b*c)-(a)^2+(0*(a+b+c+d+e))}
spp90<-function(a,b,c,d,e) {(b*c)-(b)^2+(0*(a+b+c+d+e))}

spp91<-function(a,b,c,d,e) {(d*a)-(a-b)^2 +(0*(a+b+c+d+e))}
spp92<-function(a,b,c,d,e) {(d*a)-(a-c)^2+(0*(a+b+c+d+e))}
spp93<-function(a,b,c,d,e) {(d*a)-(a-d)^2+(0*(a+b+c+d+e))}
spp94<-function(a,b,c,d,e) {(d*a)-(a-e)^2+(0*(a+b+c+d+e))}
spp95<-function(a,b,c,d,e) {(d*a)-(a+b/d)^2+(0*(a+b+c+d+e))}
spp96<-function(a,b,c,d,e) {(d*a)-(a+c/a)^2+(0*(a+b+c+d+e))}
spp97<-function(a,b,c,d,e) {(d*a)-(a+d/a)^2+(0*(a+b+c+d+e))}
spp98<-function(a,b,c,d,e) {(d*a)-(a/c+e)^2+(0*(a+b+c+d+e))}
spp99<-function(a,b,c,d,e) {(d*a)-(b)^2+(0*(a+b+c+d+e))}
spp100<-function(a,b,c,d,e) {(d*a)-(a)^2+(0*(a+b+c+d+e))}

spp101<-function(a,b,c,d,e) {(c*d*e)-(b-a)^2 +(0*(a+b+c+d+e))}
spp102<-function(a,b,c,d,e) {(c*d*e)-(b-c)^2+(0*(a+b+c+d+e))}
spp103<-function(a,b,c,d,e) {(c*d*e)-(b-d)^2+(0*(a+b+c+d+e))}
spp104<-function(a,b,c,d,e) {(c*d*e)-(b-e)^2+(0*(a+b+c+d+e))}
spp105<-function(a,b,c,d,e) {(c*d*e)-(b)^2+(0*(a+b+c+d+e))}
spp106<-function(a,b,c,d,e) {(c*d*e)-(b+c)^2+(0*(a+b+c+d+e))}
spp107<-function(a,b,c,d,e) {(c*d*e)-(b+d)^2+(0*(a+b+c+d+e))}
spp108<-function(a,b,c,d,e) {(c*d*e)-(b+e)^2+(0*(a+b+c+d+e))}
spp109<-function(a,b,c,d,e) {(c*d*e)-(b)^2+(0*(a+b+c+d+e))}
spp110<-function(a,b,c,d,e) {(c*d*e)-(a)^2+(0*(a+b+c+d+e))}

spp111<-function(a,b,c,d,e) {(d*e)-(b-a)^2 +(0*(a+b+c+d+e))}
spp112<-function(a,b,c,d,e) {(d*e)-(b-c)^2+(0*(a+b+c+d+e))}
spp113<-function(a,b,c,d,e) {(d*e)-(b-d)^2+(0*(a+b+c+d+e))}
spp114<-function(a,b,c,d,e) {(d*e)-(b-e)^2+(0*(a+b+c+d+e))}
spp115<-function(a,b,c,d,e) {(d*e)-(b+b)^2+(0*(a+b+c+d+e))}
spp116<-function(a,b,c,d,e) {(d*e)-(b+c)^2+(0*(a+b+c+d+e))}
spp117<-function(a,b,c,d,e) {(d*e)-(b/d+d)^2+(0*(a+b+c+d+e))}
spp118<-function(a,b,c,d,e) {(d*e)-(b/e+e)^2+(0*(a+b+c+d+e))}
spp119<-function(a,b,c,d,e) {(d*e)-(2*b)^2+(0*(a+b+c+d+e))}
spp120<-function(a,b,c,d,e) {(d*e)-(2*a)^2+(0*(a+b+c+d+e))}

spp121<-function(a,b,c,d,e) {(a*b)-(b-a)^2 +(0*(a+b+c+d+e))}
spp122<-function(a,b,c,d,e) {(a*b)-(b-c)^2+(0*(a+b+c+d+e))}
spp123<-function(a,b,c,d,e) {(a*b)-(b-d)^2+(0*(a+b+c+d+e))}
spp124<-function(a,b,c,d,e) {(a*b)-(b-e)^2+(0*(a+b+c+d+e))}
spp125<-function(a,b,c,d,e) {(a*b)-(b+b/c)^2+(0*(a+b+c+d+e))}
spp126<-function(a,b,c,d,e) {(a*b)-(b+c)^2+(0*(a+b+c+d+e))}
spp127<-function(a,b,c,d,e) {(a*b)-(b+d)^2+(0*(a+b+c+d+e))}
spp128<-function(a,b,c,d,e) {(a*b)-(b+e)^2+(0*(a+b+c+d+e))}
spp129<-function(a,b,c,d,e) {(a*b)-(c)^2+(0*(a+b+c+d+e))}
spp130<-function(a,b,c,d,e) {(a*b)-(d)^2+(0*(a+b+c+d+e))}

spp131<-function(a,b,c,d,e) {(c*e)-(b-a)^2 +(0*(a+b+c+d+e))}
spp132<-function(a,b,c,d,e) {(c*e)-(b-c)^2+(0*(a+b+c+d+e))}
spp133<-function(a,b,c,d,e) {(c*e)-(b-d)^2+(0*(a+b+c+d+e))}
spp134<-function(a,b,c,d,e) {(c*e)-(b-e)^2+(0*(a+b+c+d+e))}
spp135<-function(a,b,c,d,e) {(c*e)-(b+b)^2+(0*(a+b+c+d+e))}
spp136<-function(a,b,c,d,e) {(c*e)-(b+c/a)^2+(0*(a+b+c+d+e))}
spp137<-function(a,b,c,d,e) {(c*e)-(b+d)^2+(0*(a+b+c+d+e))}
spp138<-function(a,b,c,d,e) {(c*e)-(b/a+e/d)^2+(0*(a+b+c+d+e))}
spp139<-function(a,b,c,d,e) {(c*e)-(c)^2+(0*(a+b+c+d+e))}
spp140<-function(a,b,c,d,e) {(c*e)-(d)^2+(0*(a+b+c+d+e))}

spp141<-function(a,b,c,d,e) {(c*d)-(b-a)^2 +(0*(a+b+c+d+e))}
spp142<-function(a,b,c,d,e) {(c*d)-(b-c)^2+(0*(a+b+c+d+e))}
spp143<-function(a,b,c,d,e) {(c*d)-(b-d)^2+(0*(a+b+c+d+e))}
spp144<-function(a,b,c,d,e) {(c*d)-(b-e/c)^2+(0*(a+b+c+d+e))}
spp145<-function(a,b,c,d,e) {(c*d)-(b+b)^2+(0*(a+b+c+d+e))}
spp146<-function(a,b,c,d,e) {(c*d)-(b+c/a)^2+(0*(a+b+c+d+e))}
spp147<-function(a,b,c,d,e) {(c*d)-(b+d/c)^2+(0*(a+b+c+d+e))}
spp148<-function(a,b,c,d,e) {(c*d)-(b+e/b)^2+(0*(a+b+c+d+e))}
spp149<-function(a,b,c,d,e) {(c*d)-(c)^2+(0*(a+b+c+d+e))}
spp150<-function(a,b,c,d,e) {(c*d)-(d)^2+(0*(a+b+c+d+e))}

spp151<-function(a,b,c,d,e) {(b*e)-(b-a)^2 +(0*(a+b+c+d+e))}
spp152<-function(a,b,c,d,e) {(b*e)-(b-c)^2+(0*(a+b+c+d+e))}
spp153<-function(a,b,c,d,e) {(b*e)-(b-d)^2+(0*(a+b+c+d+e))}
spp154<-function(a,b,c,d,e) {(b*e)-(b-e)^2+(0*(a+b+c+d+e))}
spp155<-function(a,b,c,d,e) {(b*e)-(b+b)^2+(0*(a+b+c+d+e))}
spp156<-function(a,b,c,d,e) {(b*e)-(b+c)^2+(0*(a+b+c+d+e))}
spp157<-function(a,b,c,d,e) {(b*e)-(b+d)^2+(0*(a+b+c+d+e))}
spp158<-function(a,b,c,d,e) {(b*e)-(b/c+e/c)^2+(0*(a+b+c+d+e))}
spp159<-function(a,b,c,d,e) {(b*e)-(c)^2+(0*(a+b+c+d+e))}
spp160<-function(a,b,c,d,e) {(b*e)-(d)^2+(0*(a+b+c+d+e))}

spp161<-function(a,b,c,d,e) {(a*e)-(b-a)^2 +(0*(a+b+c+d+e))}
spp162<-function(a,b,c,d,e) {(a*e)-(b-c)^2+(0*(a+b+c+d+e))}
spp163<-function(a,b,c,d,e) {(a*e)-(b-d)^2+(0*(a+b+c+d+e))}
spp164<-function(a,b,c,d,e) {(a*e)-(b-e)^2+(0*(a+b+c+d+e))}
spp165<-function(a,b,c,d,e) {(a*e)-(b+b)^2+(0*(a+b+c+d+e))}
spp166<-function(a,b,c,d,e) {(a*e)-(b+c)^2+(0*(a+b+c+d+e))}
spp167<-function(a,b,c,d,e) {(a*e)-(b+d)^2+(0*(a+b+c+d+e))}
spp168<-function(a,b,c,d,e) {(a*e)-(b/a+e/a)^2+(0*(a+b+c+d+e))}
spp169<-function(a,b,c,d,e) {(a*e)-(c)^2+(0*(a+b+c+d+e))}
spp170<-function(a,b,c,d,e) {(a*e)-(d)^2+(0*(a+b+c+d+e))}

spp171<-function(a,b,c,d,e) {(a*c)-(b-a)^2 +(0*(a+b+c+d+e))}
spp172<-function(a,b,c,d,e) {(a*c)-(b-c)^2+(0*(a+b+c+d+e))}
spp173<-function(a,b,c,d,e) {(a*c)-(b-d)^2+(0*(a+b+c+d+e))}
spp174<-function(a,b,c,d,e) {(a*c)-(b-e/a)^2+(0*(a+b+c+d+e))}
spp175<-function(a,b,c,d,e) {(a*c)-(b+b)^2+(0*(a+b+c+d+e))}
spp176<-function(a,b,c,d,e) {(a*c)-(b/a+c)^2+(0*(a+b+c+d+e))}
spp177<-function(a,b,c,d,e) {(a*c)-(b/d+d)^2+(0*(a+b+c+d+e))}
spp178<-function(a,b,c,d,e) {(a*c)-(b/d+e)^2+(0*(a+b+c+d+e))}
spp179<-function(a,b,c,d,e) {(a*c)-(c)^2+(0*(a+b+c+d+e))}
spp180<-function(a,b,c,d,e) {(a*c)-(d)^2+(0*(a+b+c+d+e))}

spp181<-function(a,b,c,d,e) {(b*c)-(b-a)^2 +(0*(a+b+c+d+e))}
spp182<-function(a,b,c,d,e) {(b*c)-(b-c)^2+(0*(a+b+c+d+e))}
spp183<-function(a,b,c,d,e) {(b*c)-(b-d)^2+(0*(a+b+c+d+e))}
spp184<-function(a,b,c,d,e) {(b*c)-(b-e)^2+(0*(a+b+c+d+e))}
spp185<-function(a,b,c,d,e) {(b*c)-(b/c+b/c)^2+(0*(a+b+c+d+e))}
spp186<-function(a,b,c,d,e) {(b*c)-(b/a+c/a)^2+(0*(a+b+c+d+e))}
spp187<-function(a,b,c,d,e) {(b*c)-(b+d/a)^2+(0*(a+b+c+d+e))}
spp188<-function(a,b,c,d,e) {(b*c)-(b/c+e/d)^2+(0*(a+b+c+d+e))}
spp189<-function(a,b,c,d,e) {(b*c)-(c)^2+(0*(a+b+c+d+e))}
spp190<-function(a,b,c,d,e) {(b*c)-(d)^2+(0*(a+b+c+d+e))}

spp191<-function(a,b,c,d,e) {(d*a)-(b-a)^2 +(0*(a+b+c+d+e))}
spp192<-function(a,b,c,d,e) {(d*a)-(b-c)^2+(0*(a+b+c+d+e))}
spp193<-function(a,b,c,d,e) {(d*a)-(b-d)^2+(0*(a+b+c+d+e))}
spp194<-function(a,b,c,d,e) {(d*a)-(b/a-e)^2+(0*(a+b+c+d+e))}
spp195<-function(a,b,c,d,e) {(d*a)-(b/a+b)^2+(0*(a+b+c+d+e))}
spp196<-function(a,b,c,d,e) {(d*a)-(b/a+c)^2+(0*(a+b+c+d+e))}
spp197<-function(a,b,c,d,e) {(d*a)-(b/a+d)^2+(0*(a+b+c+d+e))}
spp198<-function(a,b,c,d,e) {(d*a)-(b/a+e)^2+(0*(a+b+c+d+e))}
spp199<-function(a,b,c,d,e) {(d*a)-(c)^2+(0*(a+b+c+d+e))}
spp200<-function(a,b,c,d,e) {(d*a)-(d)^2+(0*(a+b+c+d+e))}  #################################

spp201<-function(a,b,c,d,e) {(c*d*e)-(c-b)^2 +(0*(a+b+c+d+e))}
spp202<-function(a,b,c,d,e) {(c*d*e)-(c-a)^2+(0*(a+b+c+d+e))}
spp203<-function(a,b,c,d,e) {(c*d*e)-(c-d)^2+(0*(a+b+c+d+e))}
spp204<-function(a,b,c,d,e) {(c*d*e)-(c-e)^2+(0*(a+b+c+d+e))}
spp205<-function(a,b,c,d,e) {(c*d*e)-(c+b)^2+(0*(a+b+c+d+e))}
spp206<-function(a,b,c,d,e) {(c*d*e)-(c+c)^2+(0*(a+b+c+d+e))}
spp207<-function(a,b,c,d,e) {(c*d*e)-(c+d)^2+(0*(a+b+c+d+e))}
spp208<-function(a,b,c,d,e) {(c*d*e)-(c+e)^2+(0*(a+b+c+d+e))}
spp209<-function(a,b,c,d,e) {(c*d*e)-(0)^2+(0*(a+b+c+d+e))}
spp210<-function(a,b,c,d,e) {(c*d*e)-(50)^2+(0*(a+b+c+d+e))}

spp211<-function(a,b,c,d,e) {(d*e)-(c-b)^2 +(0*(a+b+c+d+e))}
spp212<-function(a,b,c,d,e) {(d*e)-(c-a)^2+(0*(a+b+c+d+e))}
spp213<-function(a,b,c,d,e) {(d*e)-(c-d)^2+(0*(a+b+c+d+e))}
spp214<-function(a,b,c,d,e) {(d*e)-(c-e)^2+(0*(a+b+c+d+e))}
spp215<-function(a,b,c,d,e) {(d*e)-(c/b)^2+(0*(a+b+c+d+e))}
spp216<-function(a,b,c,d,e) {(d*e)-(c/c)^2+(0*(a+b+c+d+e))}
spp217<-function(a,b,c,d,e) {(d*e)-(c/d)^2+(0*(a+b+c+d+e))}
spp218<-function(a,b,c,d,e) {(d*e)-(c/e)^2+(0*(a+b+c+d+e))}
spp219<-function(a,b,c,d,e) {(d*e)-(0)^2+(0*(a+b+c+d+e))}
spp220<-function(a,b,c,d,e) {(d*e)-(5)^2+(0*(a+b+c+d+e))}

spp221<-function(a,b,c,d,e) {(a*b)-(c-b)^2 +(0*(a+b+c+d+e))}
spp222<-function(a,b,c,d,e) {(a*b)-(c-a)^2+(0*(a+b+c+d+e))}
spp223<-function(a,b,c,d,e) {(a*b)-(c-d)^2+(0*(a+b+c+d+e))}
spp224<-function(a,b,c,d,e) {(a*b)-(c-e)^2+(0*(a+b+c+d+e))}
spp225<-function(a,b,c,d,e) {(a*b)-(a/b)^2+(0*(a+b+c+d+e))}
spp226<-function(a,b,c,d,e) {(a*b)-(a/c)^2+(0*(a+b+c+d+e))}
spp227<-function(a,b,c,d,e) {(a*b)-(a/d)^2+(0*(a+b+c+d+e))}
spp228<-function(a,b,c,d,e) {(a*b)-(a/e)^2+(0*(a+b+c+d+e))}
spp229<-function(a,b,c,d,e) {(a*b)-(5)^2+(0*(a+b+c+d+e))}
spp230<-function(a,b,c,d,e) {(a*b)-(50)^2+(0*(a+b+c+d+e))}

spp231<-function(a,b,c,d,e) {(c*e)-(c-b)^2 +(0*(a+b+c+d+e))}
spp232<-function(a,b,c,d,e) {(c*e)-(c-a)^2+(0*(a+b+c+d+e))}
spp233<-function(a,b,c,d,e) {(c*e)-(c-d)^2+(0*(a+b+c+d+e))}
spp234<-function(a,b,c,d,e) {(c*e)-(c-e)^2+(0*(a+b+c+d+e))}
spp235<-function(a,b,c,d,e) {(c*e)-(a/b)^2+(0*(a+b+c+d+e))}
spp236<-function(a,b,c,d,e) {(c*e)-(a/c)^2+(0*(a+b+c+d+e))}
spp237<-function(a,b,c,d,e) {(c*e)-(a/d)^2+(0*(a+b+c+d+e))}
spp238<-function(a,b,c,d,e) {(c*e)-(a/e)^2+(0*(a+b+c+d+e))}
spp239<-function(a,b,c,d,e) {(c*e)-(0)^2+(0*(a+b+c+d+e))}
spp240<-function(a,b,c,d,e) {(c*e)-(5)^2+(0*(a+b+c+d+e))}

spp241<-function(a,b,c,d,e) {(c*d)-(c-b)^2 +(0*(a+b+c+d+e))}
spp242<-function(a,b,c,d,e) {(c*d)-(c-a)^2+(0*(a+b+c+d+e))}
spp243<-function(a,b,c,d,e) {(c*d)-(c-d)^2+(0*(a+b+c+d+e))}
spp244<-function(a,b,c,d,e) {(c*d)-(c-e)^2+(0*(a+b+c+d+e))}
spp245<-function(a,b,c,d,e) {(c*d)-(a/b)^2+(0*(a+b+c+d+e))}
spp246<-function(a,b,c,d,e) {(c*d)-(a/c)^2+(0*(a+b+c+d+e))}
spp247<-function(a,b,c,d,e) {(c*d)-(a/d)^2+(0*(a+b+c+d+e))}
spp248<-function(a,b,c,d,e) {(c*d)-(a/e)^2+(0*(a+b+c+d+e))}
spp249<-function(a,b,c,d,e) {(c*d)-(0)^2+(0*(a+b+c+d+e))}
spp250<-function(a,b,c,d,e) {(c*d)-(5)^2+(0*(a+b+c+d+e))}

spp251<-function(a,b,c,d,e) {(b*e)-(c-b)^2 +(0*(a+b+c+d+e))}
spp252<-function(a,b,c,d,e) {(b*e)-(c-a)^2+(0*(a+b+c+d+e))}
spp253<-function(a,b,c,d,e) {(b*e)-(c-d)^2+(0*(a+b+c+d+e))}
spp254<-function(a,b,c,d,e) {(b*e)-(c-e)^2+(0*(a+b+c+d+e))}
spp255<-function(a,b,c,d,e) {(b*e)-(a/b)^2+(0*(a+b+c+d+e))}
spp256<-function(a,b,c,d,e) {(b*e)-(a/c)^2+(0*(a+b+c+d+e))}
spp257<-function(a,b,c,d,e) {(b*e)-(a/d)^2+(0*(a+b+c+d+e))}
spp258<-function(a,b,c,d,e) {(b*e)-(a/e)^2+(0*(a+b+c+d+e))}
spp259<-function(a,b,c,d,e) {(b*e)-(0)^2+(0*(a+b+c+d+e))}
spp260<-function(a,b,c,d,e) {(b*e)-(50)^2+(0*(a+b+c+d+e))}

spp261<-function(a,b,c,d,e) {(a*e)-(c-b)^2 +(0*(a+b+c+d+e))}
spp262<-function(a,b,c,d,e) {(a*e)-(c-a)^2+(0*(a+b+c+d+e))}
spp263<-function(a,b,c,d,e) {(a*e)-(c-d)^2+(0*(a+b+c+d+e))}
spp264<-function(a,b,c,d,e) {(a*e)-(c-e)^2+(0*(a+b+c+d+e))}
spp265<-function(a,b,c,d,e) {(a*e)-(a/b)^2+(0*(a+b+c+d+e))}
spp266<-function(a,b,c,d,e) {(a*e)-(a/c)^2+(0*(a+b+c+d+e))}
spp267<-function(a,b,c,d,e) {(a*e)-(a/d)^2+(0*(a+b+c+d+e))}
spp268<-function(a,b,c,d,e) {(a*e)-(a/e)^2+(0*(a+b+c+d+e))}
spp269<-function(a,b,c,d,e) {(a*e)-(0)^2+(0*(a+b+c+d+e))}
spp270<-function(a,b,c,d,e) {(a*e)-(5)^2+(0*(a+b+c+d+e))}

spp271<-function(a,b,c,d,e) {(a*c)-(c-b)^2 +(0*(a+b+c+d+e))}
spp272<-function(a,b,c,d,e) {(a*c)-(c-a)^2+(0*(a+b+c+d+e))}
spp273<-function(a,b,c,d,e) {(a*c)-(c-d)^2+(0*(a+b+c+d+e))}
spp274<-function(a,b,c,d,e) {(a*c)-(c-e)^2+(0*(a+b+c+d+e))}
spp275<-function(a,b,c,d,e) {(a*c)-(a/b)^2+(0*(a+b+c+d+e))}
spp276<-function(a,b,c,d,e) {(a*c)-(a/c)^2+(0*(a+b+c+d+e))}
spp277<-function(a,b,c,d,e) {(a*c)-(a/d)^2+(0*(a+b+c+d+e))}
spp278<-function(a,b,c,d,e) {(a*c)-(a/e)^2+(0*(a+b+c+d+e))}
spp279<-function(a,b,c,d,e) {(a*c)-(0)^2+(0*(a+b+c+d+e))}
spp280<-function(a,b,c,d,e) {(a*c)-(5)^2+(0*(a+b+c+d+e))}

spp281<-function(a,b,c,d,e) {(b*c)-(c-b)^2 +(0*(a+b+c+d+e))}
spp282<-function(a,b,c,d,e) {(b*c)-(c-a)^2+(0*(a+b+c+d+e))}
spp283<-function(a,b,c,d,e) {(b*c)-(c-d)^2+(0*(a+b+c+d+e))}
spp284<-function(a,b,c,d,e) {(b*c)-(c-e)^2+(0*(a+b+c+d+e))}
spp285<-function(a,b,c,d,e) {(b*c)-(a/b)^2+(0*(a+b+c+d+e))}
spp286<-function(a,b,c,d,e) {(b*c)-(a/c)^2+(0*(a+b+c+d+e))}
spp287<-function(a,b,c,d,e) {(b*c)-(a/d)^2+(0*(a+b+c+d+e))}
spp288<-function(a,b,c,d,e) {(b*c)-(a/e)^2+(0*(a+b+c+d+e))}
spp289<-function(a,b,c,d,e) {(b*c)-(0)^2+(0*(a+b+c+d+e))}
spp290<-function(a,b,c,d,e) {(b*c)-(5)^2+(0*(a+b+c+d+e))}

spp291<-function(a,b,c,d,e) {(d*a)-(c-b)^2 +(0*(a+b+c+d+e))}
spp292<-function(a,b,c,d,e) {(d*a)-(c-a)^2+(0*(a+b+c+d+e))}
spp293<-function(a,b,c,d,e) {(d*a)-(c-d)^2+(0*(a+b+c+d+e))}
spp294<-function(a,b,c,d,e) {(d*a)-(c-e)^2+(0*(a+b+c+d+e))}
spp295<-function(a,b,c,d,e) {(d*a)-(a/b)^2+(0*(a+b+c+d+e))}
spp296<-function(a,b,c,d,e) {(d*a)-(a/c)^2+(0*(a+b+c+d+e))}
spp297<-function(a,b,c,d,e) {(d*a)-(a/d)^2+(0*(a+b+c+d+e))}
spp298<-function(a,b,c,d,e) {(d*a)-(a/e)^2+(0*(a+b+c+d+e))}
spp299<-function(a,b,c,d,e) {(d*a)-(0)^2+(0*(a+b+c+d+e))}
spp300<-function(a,b,c,d,e) {(d*a)-(5)^2+(0*(a+b+c+d+e))}###################################

spp301<-function(a,b,c,d,e) {(100*a)-(d-b)^2 +(0*(a+b+c+d+e))}
spp302<-function(a,b,c,d,e) {(100*a)-(d-c)^2+(0*(a+b+c+d+e))}
spp303<-function(a,b,c,d,e) {(100*a)-(d-a)^2+(0*(a+b+c+d+e))}
spp304<-function(a,b,c,d,e) {(100*a)-(d-e)^2+(0*(a+b+c+d+e))}
spp305<-function(a,b,c,d,e) {(100*a)-(d+b)^2+(0*(a+b+c+d+e))}
spp306<-function(a,b,c,d,e) {(100*a)-(d+c)^2+(0*(a+b+c+d+e))}
spp307<-function(a,b,c,d,e) {(100*a)-(d+d)^2+(0*(a+b+c+d+e))}
spp308<-function(a,b,c,d,e) {(100*a)-(d+e)^2+(0*(a+b+c+d+e))}
spp309<-function(a,b,c,d,e) {(1000*a)-(100)^2+(0*(a+b+c+d+e))}
spp310<-function(a,b,c,d,e) {(10000*a)-(100)^2+(0*(a+b+c+d+e))}

spp311<-function(a,b,c,d,e) {(100*b)-(d-b)^2 +(0*(a+b+c+d+e))}
spp312<-function(a,b,c,d,e) {(100*b)-(d-c)^2+(0*(a+b+c+d+e))}
spp313<-function(a,b,c,d,e) {(100*b)-(d-a)^2+(0*(a+b+c+d+e))}
spp314<-function(a,b,c,d,e) {(100*b)-(d-e)^2+(0*(a+b+c+d+e))}
spp315<-function(a,b,c,d,e) {(100*b)-(d+b)^2+(0*(a+b+c+d+e))}
spp316<-function(a,b,c,d,e) {(100*b)-(d+c)^2+(0*(a+b+c+d+e))}
spp317<-function(a,b,c,d,e) {(100*b)-(d+d)^2+(0*(a+b+c+d+e))}
spp318<-function(a,b,c,d,e) {(100*b)-(d+e)^2+(0*(a+b+c+d+e))}
spp319<-function(a,b,c,d,e) {(1000*b)-(10)^2+(0*(a+b+c+d+e))}
spp320<-function(a,b,c,d,e) {(10000*b)-(100)^2+(0*(a+b+c+d+e))}

spp321<-function(a,b,c,d,e) {(100*c)-(d-b)^2 +(0*(a+b+c+d+e))}
spp322<-function(a,b,c,d,e) {(100*c)-(d-c)^2+(0*(a+b+c+d+e))}
spp323<-function(a,b,c,d,e) {(100*c)-(d-a)^2+(0*(a+b+c+d+e))}
spp324<-function(a,b,c,d,e) {(100*c)-(d-e)^2+(0*(a+b+c+d+e))}
spp325<-function(a,b,c,d,e) {(100*c)-(d+b)^2+(0*(a+b+c+d+e))}
spp326<-function(a,b,c,d,e) {(100*c)-(d+c)^2+(0*(a+b+c+d+e))}
spp327<-function(a,b,c,d,e) {(100*c)-(d+d)^2+(0*(a+b+c+d+e))}
spp328<-function(a,b,c,d,e) {(100*c)-(d+e)^2+(0*(a+b+c+d+e))}
spp329<-function(a,b,c,d,e) {(1000*c)-(10)^2+(0*(a+b+c+d+e))}
spp330<-function(a,b,c,d,e) {(1000*c)-(10)^2+(0*(a+b+c+d+e))}

spp331<-function(a,b,c,d,e) {(100*d)-(d-b)^2 +(0*(a+b+c+d+e))}
spp332<-function(a,b,c,d,e) {(100*d)-(d-c)^2+(0*(a+b+c+d+e))}
spp333<-function(a,b,c,d,e) {(100*d)-(d-a)^2+(0*(a+b+c+d+e))}
spp334<-function(a,b,c,d,e) {(100*d)-(d-e)^2+(0*(a+b+c+d+e))}
spp335<-function(a,b,c,d,e) {(100*d)-(d+b)^2+(0*(a+b+c+d+e))}
spp336<-function(a,b,c,d,e) {(100*d)-(d+c)^2+(0*(a+b+c+d+e))}
spp337<-function(a,b,c,d,e) {(100*d)-(d+d)^2+(0*(a+b+c+d+e))}
spp338<-function(a,b,c,d,e) {(100*d)-(d+e)^2+(0*(a+b+c+d+e))}
spp339<-function(a,b,c,d,e) {(10000*d)-(10)^2+(0*(a+b+c+d+e))}
spp340<-function(a,b,c,d,e) {(10000*d)-(100)^2+(0*(a+b+c+d+e))}

spp341<-function(a,b,c,d,e) {(300*d)-(d-b)^2 +(0*(a+b+c+d+e))}
spp342<-function(a,b,c,d,e) {(300*d)-(d-c)^2+(0*(a+b+c+d+e))}
spp343<-function(a,b,c,d,e) {(300*d)-(d-a)^2+(0*(a+b+c+d+e))}
spp344<-function(a,b,c,d,e) {(300*d)-(d-e)^2+(0*(a+b+c+d+e))}
spp345<-function(a,b,c,d,e) {(300*d)-(d+b)^2+(0*(a+b+c+d+e))}
spp346<-function(a,b,c,d,e) {(300*d)-(d+c)^2+(0*(a+b+c+d+e))}
spp347<-function(a,b,c,d,e) {(300*d)-(d+d)^2+(0*(a+b+c+d+e))}
spp348<-function(a,b,c,d,e) {(300*d)-(d+e)^2+(0*(a+b+c+d+e))}
spp349<-function(a,b,c,d,e) {(3000*d)-(100)^2+(0*(a+b+c+d+e))}
spp350<-function(a,b,c,d,e) {(30000*d)-(100)^2+(0*(a+b+c+d+e))}

spp351<-function(a,b,c,d,e) {(100*e)-(d-b)^2 +(0*(a+b+c+d+e))}
spp352<-function(a,b,c,d,e) {(100*e)-(d-c)^2+(0*(a+b+c+d+e))}
spp353<-function(a,b,c,d,e) {(100*e)-(d-a)^2+(0*(a+b+c+d+e))}
spp354<-function(a,b,c,d,e) {(100*e)-(d-e)^2+(0*(a+b+c+d+e))}
spp355<-function(a,b,c,d,e) {(100*e)-(d+b)^2+(0*(a+b+c+d+e))}
spp356<-function(a,b,c,d,e) {(100*e)-(d+c)^2+(0*(a+b+c+d+e))}
spp357<-function(a,b,c,d,e) {(100*e)-(d+d)^2+(0*(a+b+c+d+e))}
spp358<-function(a,b,c,d,e) {(100*e)-(d+e)^2+(0*(a+b+c+d+e))}
spp359<-function(a,b,c,d,e) {(100*e)-(10)^2+(0*(a+b+c+d+e))}
spp360<-function(a,b,c,d,e) {(1000*e)-(100)^2+(0*(a+b+c+d+e))}

spp361<-function(a,b,c,d,e) {(1000*e)-(d-b)^2 +(0*(a+b+c+d+e))}
spp362<-function(a,b,c,d,e) {(1000*e)-(d-c)^2+(0*(a+b+c+d+e))}
spp363<-function(a,b,c,d,e) {(1000*e)-(d-a)^2+(0*(a+b+c+d+e))}
spp364<-function(a,b,c,d,e) {(1000*e)-(d-e)^2+(0*(a+b+c+d+e))}
spp365<-function(a,b,c,d,e) {(1000*e)-(d+b)^2+(0*(a+b+c+d+e))}
spp366<-function(a,b,c,d,e) {(1000*e)-(d+c)^2+(0*(a+b+c+d+e))}
spp367<-function(a,b,c,d,e) {(1000*e)-(d+d)^2+(0*(a+b+c+d+e))}
spp368<-function(a,b,c,d,e) {(1000*e)-(d+e)^2+(0*(a+b+c+d+e))}
spp369<-function(a,b,c,d,e) {(1000*e)-(100)^2+(0*(a+b+c+d+e))}
spp370<-function(a,b,c,d,e) {(1000*e)-(200)^2+(0*(a+b+c+d+e))}

spp371<-function(a,b,c,d,e) {(1000*d)-(d-b)^2 +(0*(a+b+c+d+e))}
spp372<-function(a,b,c,d,e) {(1000*d)-(d-c)^2+(0*(a+b+c+d+e))}
spp373<-function(a,b,c,d,e) {(1000*d)-(d-a)^2+(0*(a+b+c+d+e))}
spp374<-function(a,b,c,d,e) {(1000*d)-(d-e)^2+(0*(a+b+c+d+e))}
spp375<-function(a,b,c,d,e) {(1000*d)-(d+b)^2+(0*(a+b+c+d+e))}
spp376<-function(a,b,c,d,e) {(1000*d)-(d+c)^2+(0*(a+b+c+d+e))}
spp377<-function(a,b,c,d,e) {(1000*d)-(d+d)^2+(0*(a+b+c+d+e))}
spp378<-function(a,b,c,d,e) {(1000*d)-(d+e)^2+(0*(a+b+c+d+e))}
spp379<-function(a,b,c,d,e) {(1000*d)-(100)^2+(0*(a+b+c+d+e))}

spp380<-function(a,b,c,d,e) {(1000*d)-(10)^2+(0*(a+b+c+d+e))}
spp381<-function(a,b,c,d,e) {(1000*c)-(d-b)^2 +(0*(a+b+c+d+e))}
spp382<-function(a,b,c,d,e) {(1000*c)-(d-c)^2+(0*(a+b+c+d+e))}
spp383<-function(a,b,c,d,e) {(1000*c)-(d-a)^2+(0*(a+b+c+d+e))}
spp384<-function(a,b,c,d,e) {(1000*c)-(d-e)^2+(0*(a+b+c+d+e))}
spp385<-function(a,b,c,d,e) {(1000*c)-(d+b)^2+(0*(a+b+c+d+e))}
spp386<-function(a,b,c,d,e) {(1000*c)-(d+c)^2+(0*(a+b+c+d+e))}
spp387<-function(a,b,c,d,e) {(1000*c)-(d+d)^2+(0*(a+b+c+d+e))}
spp388<-function(a,b,c,d,e) {(1000*c)-(d+e)^2+(0*(a+b+c+d+e))}
spp389<-function(a,b,c,d,e) {(1000*c)-(100)^2+(0*(a+b+c+d+e))}
spp390<-function(a,b,c,d,e) {(1000*c)-(200/a)^2+(0*(a+b+c+d+e))}

spp391<-function(a,b,c,d,e) {(1000*b)-(d-b)^2 +(0*(a+b+c+d+e))}
spp392<-function(a,b,c,d,e) {(1000*b)-(d-c)^2+(0*(a+b+c+d+e))}
spp393<-function(a,b,c,d,e) {(1000*b)-(d-a)^2+(0*(a+b+c+d+e))}
spp394<-function(a,b,c,d,e) {(1000*b)-(d-e)^2+(0*(a+b+c+d+e))}
spp395<-function(a,b,c,d,e) {(1000*b)-(d+b)^2+(0*(a+b+c+d+e))}
spp396<-function(a,b,c,d,e) {(1000*b)-(d+c)^2+(0*(a+b+c+d+e))}
spp397<-function(a,b,c,d,e) {(1000*b)-(d+d)^2+(0*(a+b+c+d+e))}
spp398<-function(a,b,c,d,e) {(1000*b)-(d+e)^2+(0*(a+b+c+d+e))}
spp399<-function(a,b,c,d,e) {(1000*b)-(100)^2+(0*(a+b+c+d+e))}
spp400<-function(a,b,c,d,e) {(1000*b)-(200)^2+(0*(a+b+c+d+e))}

spp401<-function(a,b,c,d,e) {(c*d*e)-(d-b)^2 +(0*(a+b+c+d+e))}
spp402<-function(a,b,c,d,e) {(c*d*e)-(d-c)^2+(0*(a+b+c+d+e))}
spp403<-function(a,b,c,d,e) {(c*d*e)-(d-a)^2+(0*(a+b+c+d+e))}
spp404<-function(a,b,c,d,e) {(c*d*e)-(d-e)^2+(0*(a+b+c+d+e))}
spp405<-function(a,b,c,d,e) {(c*d*e)-(d+b)^2+(0*(a+b+c+d+e))}
spp406<-function(a,b,c,d,e) {(c*d*e)-(d+c)^2+(0*(a+b+c+d+e))}
spp407<-function(a,b,c,d,e) {(c*d*e)-(d+d)^2+(0*(a+b+c+d+e))}
spp408<-function(a,b,c,d,e) {(c*d*e)-(d+e)^2+(0*(a+b+c+d+e))}
spp409<-function(a,b,c,d,e) {(c*d*e)-(100)^2+(0*(a+b+c+d+e))}
spp410<-function(a,b,c,d,e) {(c*d*e)-(200/a)^2+(0*(a+b+c+d+e))}

spp411<-function(a,b,c,d,e) {(d*e)-(d-b)^2 +(0*(a+b+c+d+e))}
spp412<-function(a,b,c,d,e) {(d*e)-(d-c)^2+(0*(a+b+c+d+e))}
spp413<-function(a,b,c,d,e) {(d*e)-(d-a)^2+(0*(a+b+c+d+e))}
spp414<-function(a,b,c,d,e) {(d*e)-(d-e)^2+(0*(a+b+c+d+e))}
spp415<-function(a,b,c,d,e) {(d*e)-(d/c+b/c)^2+(0*(a+b+c+d+e))}
spp416<-function(a,b,c,d,e) {(d*e)-(d+c)^2+(0*(a+b+c+d+e))}
spp417<-function(a,b,c,d,e) {(d*e)-(d+d)^2+(0*(a+b+c+d+e))}
spp418<-function(a,b,c,d,e) {(d*e)-(d/e+e/d)^2+(0*(a+b+c+d+e))}
spp419<-function(a,b,c,d,e) {(d*e)-(10)^2+(0*(a+b+c+d+e))}
spp420<-function(a,b,c,d,e) {(d*e)-(20)^2+(0*(a+b+c+d+e))}

spp421<-function(a,b,c,d,e) {(a*b)-(d-b)^2 +(0*(a+b+c+d+e))}
spp422<-function(a,b,c,d,e) {(a*b)-(d-c)^2+(0*(a+b+c+d+e))}
spp423<-function(a,b,c,d,e) {(a*b)-(d-a)^2+(0*(a+b+c+d+e))}
spp424<-function(a,b,c,d,e) {(a*b)-(d-e)^2+(0*(a+b+c+d+e))}
spp425<-function(a,b,c,d,e) {(a*b)-(d+b)^2+(0*(a+b+c+d+e))}
spp426<-function(a,b,c,d,e) {(a*b)-(d+c)^2+(0*(a+b+c+d+e))}
spp427<-function(a,b,c,d,e) {(a*b)-(d+d)^2+(0*(a+b+c+d+e))}
spp428<-function(a,b,c,d,e) {(a*b)-(d+e)^2+(0*(a+b+c+d+e))}
spp429<-function(a,b,c,d,e) {(a*b)-(10)^2+(0*(a+b+c+d+e))}
spp430<-function(a,b,c,d,e) {(a*b)-(50)^2+(0*(a+b+c+d+e))}

spp431<-function(a,b,c,d,e) {(c*e)-(d-b)^2 +(0*(a+b+c+d+e))}
spp432<-function(a,b,c,d,e) {(c*e)-(d-c)^2+(0*(a+b+c+d+e))}
spp433<-function(a,b,c,d,e) {(c*e)-(d-a)^2+(0*(a+b+c+d+e))}
spp434<-function(a,b,c,d,e) {(c*e)-(d-e)^2+(0*(a+b+c+d+e))}
spp435<-function(a,b,c,d,e) {(c*e)-(d+b)^2+(0*(a+b+c+d+e))}
spp436<-function(a,b,c,d,e) {(c*e)-(d+c)^2+(0*(a+b+c+d+e))}
spp437<-function(a,b,c,d,e) {(c*e)-(d+d)^2+(0*(a+b+c+d+e))}
spp438<-function(a,b,c,d,e) {(c*e)-(d/c+e/c)^2+(0*(a+b+c+d+e))}
spp439<-function(a,b,c,d,e) {(c*e)-(10)^2+(0*(a+b+c+d+e))}
spp440<-function(a,b,c,d,e) {(c*e)-(200/d)^2+(0*(a+b+c+d+e))}

spp441<-function(a,b,c,d,e) {(c*d)-(d-b)^2 +(0*(a+b+c+d+e))}
spp442<-function(a,b,c,d,e) {(c*d)-(d-c)^2+(0*(a+b+c+d+e))}
spp443<-function(a,b,c,d,e) {(c*d)-(d-a)^2+(0*(a+b+c+d+e))}
spp444<-function(a,b,c,d,e) {(c*d)-(d-e)^2+(0*(a+b+c+d+e))}
spp445<-function(a,b,c,d,e) {(c*d)-(d/a+b/a)^2+(0*(a+b+c+d+e))}
spp446<-function(a,b,c,d,e) {(c*d)-(d/a+c/a)^2+(0*(a+b+c+d+e))}
spp447<-function(a,b,c,d,e) {(c*d)-(d/c+d/a)^2+(0*(a+b+c+d+e))}
spp448<-function(a,b,c,d,e) {(c*d)-(d/a+e/a)^2+(0*(a+b+c+d+e))}
spp449<-function(a,b,c,d,e) {(c*d)-(10/a)^2+(0*(a+b+c+d+e))}
spp450<-function(a,b,c,d,e) {(c*d)-(100/a)^2+(0*(a+b+c+d+e))}

spp451<-function(a,b,c,d,e) {(b*e)-(d-b)^2 +(0*(a+b+c+d+e))}
spp452<-function(a,b,c,d,e) {(b*e)-(d-c)^2+(0*(a+b+c+d+e))}
spp453<-function(a,b,c,d,e) {(b*e)-(d-a)^2+(0*(a+b+c+d+e))}
spp454<-function(a,b,c,d,e) {(b*e)-(d-e)^2+(0*(a+b+c+d+e))}
spp455<-function(a,b,c,d,e) {(b*e)-(d+b)^2+(0*(a+b+c+d+e))}
spp456<-function(a,b,c,d,e) {(b*e)-(d+c)^2+(0*(a+b+c+d+e))}
spp457<-function(a,b,c,d,e) {(b*e)-(d+d)^2+(0*(a+b+c+d+e))}
spp458<-function(a,b,c,d,e) {(b*e)-(d+e)^2+(0*(a+b+c+d+e))}
spp459<-function(a,b,c,d,e) {(b*e)-(10)^2+(0*(a+b+c+d+e))}
spp460<-function(a,b,c,d,e) {(b*e)-(100/a)^2+(0*(a+b+c+d+e))}

spp461<-function(a,b,c,d,e) {(a*e)-(d-b)^2 +(0*(a+b+c+d+e))}
spp462<-function(a,b,c,d,e) {(a*e)-(d-c)^2+(0*(a+b+c+d+e))}
spp463<-function(a,b,c,d,e) {(a*e)-(d-a)^2+(0*(a+b+c+d+e))}
spp464<-function(a,b,c,d,e) {(a*e)-(d-e)^2+(0*(a+b+c+d+e))}
spp465<-function(a,b,c,d,e) {(a*e)-(d/c+b)^2+(0*(a+b+c+d+e))}
spp466<-function(a,b,c,d,e) {(a*e)-(d+c)^2+(0*(a+b+c+d+e))}
spp467<-function(a,b,c,d,e) {(a*e)-(d+d)^2+(0*(a+b+c+d+e))}
spp468<-function(a,b,c,d,e) {(a*e)-(d+e)^2+(0*(a+b+c+d+e))}
spp469<-function(a,b,c,d,e) {(a*e)-(10)^2+(0*(a+b+c+d+e))}
spp470<-function(a,b,c,d,e) {(a*e)-(100/a)^2+(0*(a+b+c+d+e))}

spp471<-function(a,b,c,d,e) {(a*c)-(d-b)^2 +(0*(a+b+c+d+e))}
spp472<-function(a,b,c,d,e) {(a*c)-(d-c)^2+(0*(a+b+c+d+e))}
spp473<-function(a,b,c,d,e) {(a*c)-(d-a)^2+(0*(a+b+c+d+e))}
spp474<-function(a,b,c,d,e) {(a*c)-(d-e)^2+(0*(a+b+c+d+e))}
spp475<-function(a,b,c,d,e) {(a*c)-(d/e+b/e)^2+(0*(a+b+c+d+e))}
spp476<-function(a,b,c,d,e) {(a*c)-(d+c)^2+(0*(a+b+c+d+e))}
spp477<-function(a,b,c,d,e) {(a*c)-(d+d)^2+(0*(a+b+c+d+e))}
spp478<-function(a,b,c,d,e) {(a*c)-(d+e)^2+(0*(a+b+c+d+e))}
spp479<-function(a,b,c,d,e) {(a*c)-(10)^2+(0*(a+b+c+d+e))}
spp480<-function(a,b,c,d,e) {(a*c)-(100/d)^2+(0*(a+b+c+d+e))}

spp481<-function(a,b,c,d,e) {(b*c)-(d-b)^2 +(0*(a+b+c+d+e))}
spp482<-function(a,b,c,d,e) {(b*c)-(d-c)^2+(0*(a+b+c+d+e))}
spp483<-function(a,b,c,d,e) {(b*c)-(d-a)^2+(0*(a+b+c+d+e))}
spp484<-function(a,b,c,d,e) {(b*c)-(d-e)^2+(0*(a+b+c+d+e))}
spp485<-function(a,b,c,d,e) {(b*c)-(d/a+b/c)^2+(0*(a+b+c+d+e))}
spp486<-function(a,b,c,d,e) {(b*c)-(d+c)^2+(0*(a+b+c+d+e))}
spp487<-function(a,b,c,d,e) {(b*c)-(d+d)^2+(0*(a+b+c+d+e))}
spp488<-function(a,b,c,d,e) {(b*c)-(d+e)^2+(0*(a+b+c+d+e))}
spp489<-function(a,b,c,d,e) {(b*c)-(10)^2+(0*(a+b+c+d+e))}
spp490<-function(a,b,c,d,e) {(b*c)-(100/a)^2+(0*(a+b+c+d+e))}

spp491<-function(a,b,c,d,e) {(d*a)-(d-b)^2 +(0*(a+b+c+d+e))}
spp492<-function(a,b,c,d,e) {(d*a)-(d-c)^2+(0*(a+b+c+d+e))}
spp493<-function(a,b,c,d,e) {(d*a)-(d-a)^2+(0*(a+b+c+d+e))}
spp494<-function(a,b,c,d,e) {(d*a)-(d-e)^2+(0*(a+b+c+d+e))}
spp495<-function(a,b,c,d,e) {(d*a)-(d/a+b/a)^2+(0*(a+b+c+d+e))}
spp496<-function(a,b,c,d,e) {(d*a)-(d+c)^2+(0*(a+b+c+d+e))}
spp497<-function(a,b,c,d,e) {(d*a)-(d+d)^2+(0*(a+b+c+d+e))}
spp498<-function(a,b,c,d,e) {(d*a)-(d+e)^2+(0*(a+b+c+d+e))}
spp499<-function(a,b,c,d,e) {(d*a)-(10)^2+(0*(a+b+c+d+e))}
spp500<-function(a,b,c,d,e) {(d*a)-(100/c)^2+(0*(a+b+c+d+e))}###################################


spp501<-function(a,b,c,d,e) {(c*d*e)-(e-b)^2 +(0*(a+b+c+d+e))}
spp502<-function(a,b,c,d,e) {(c*d*e)-(e-c)^2+(0*(a+b+c+d+e))}
spp503<-function(a,b,c,d,e) {(c*d*e)-(e-d)^2+(0*(a+b+c+d+e))}
spp504<-function(a,b,c,d,e) {(c*d*e)-(e-a)^2+(0*(a+b+c+d+e))}
spp505<-function(a,b,c,d,e) {(c*d*e)-(e+b)^2+(0*(a+b+c+d+e))}
spp506<-function(a,b,c,d,e) {(c*d*e)-(e+c)^2+(0*(a+b+c+d+e))}
spp507<-function(a,b,c,d,e) {(c*d*e)-(e+d)^2+(0*(a+b+c+d+e))}
spp508<-function(a,b,c,d,e) {(c*d*e)-(e+e)^2+(0*(a+b+c+d+e))}
spp509<-function(a,b,c,d,e) {(c*d*e)-(2)^2+(0*(a+b+c+d+e))}
spp510<-function(a,b,c,d,e) {(c*d*e)-(20)^2+(0*(a+b+c+d+e))}

spp511<-function(a,b,c,d,e) {(d*e)-(e-b)^2 +(0*(a+b+c+d+e))}
spp512<-function(a,b,c,d,e) {(d*e)-(e-c)^2+(0*(a+b+c+d+e))}
spp513<-function(a,b,c,d,e) {(d*e)-(e-d)^2+(0*(a+b+c+d+e))}
spp514<-function(a,b,c,d,e) {(d*e)-(e-a)^2+(0*(a+b+c+d+e))}
spp515<-function(a,b,c,d,e) {(d*e)-(e/a+b/a)^2+(0*(a+b+c+d+e))}
spp516<-function(a,b,c,d,e) {(d*e)-(e/a+c/a)^2+(0*(a+b+c+d+e))}
spp517<-function(a,b,c,d,e) {(d*e)-(e/a+d/a)^2+(0*(a+b+c+d+e))}
spp518<-function(a,b,c,d,e) {(d*e)-(e/a+e/a)^2+(0*(a+b+c+d+e))}
spp519<-function(a,b,c,d,e) {(d*e)-(2)^2+(0*(a+b+c+d+e))}
spp520<-function(a,b,c,d,e) {(d*e)-(20)^2+(0*(a+b+c+d+e))}

spp521<-function(a,b,c,d,e) {(a*b)-(e-b)^2 +(0*(a+b+c+d+e))}
spp522<-function(a,b,c,d,e) {(a*b)-(e-c)^2+(0*(a+b+c+d+e))}
spp523<-function(a,b,c,d,e) {(a*b)-(e-d)^2+(0*(a+b+c+d+e))}
spp524<-function(a,b,c,d,e) {(a*b)-(e-a)^2+(0*(a+b+c+d+e))}
spp525<-function(a,b,c,d,e) {(a*b)-(e+b)^2+(0*(a+b+c+d+e))}
spp526<-function(a,b,c,d,e) {(a*b)-(e+c)^2+(0*(a+b+c+d+e))}
spp527<-function(a,b,c,d,e) {(a*b)-(e+d)^2+(0*(a+b+c+d+e))}
spp528<-function(a,b,c,d,e) {(a*b)-(e+e)^2+(0*(a+b+c+d+e))}
spp529<-function(a,b,c,d,e) {(a*b)-(2)^2+(0*(a+b+c+d+e))}
spp530<-function(a,b,c,d,e) {(a*b)-(20)^2+(0*(a+b+c+d+e))}

spp531<-function(a,b,c,d,e) {(c*e)-(e-b)^2 +(0*(a+b+c+d+e))}
spp532<-function(a,b,c,d,e) {(c*e)-(e-c)^2+(0*(a+b+c+d+e))}
spp533<-function(a,b,c,d,e) {(c*e)-(e-d)^2+(0*(a+b+c+d+e))}
spp534<-function(a,b,c,d,e) {(c*e)-(e-a)^2+(0*(a+b+c+d+e))}
spp535<-function(a,b,c,d,e) {(c*e)-(e/a+b/a)^2+(0*(a+b+c+d+e))}
spp536<-function(a,b,c,d,e) {(c*e)-(e/a+c/a)^2+(0*(a+b+c+d+e))}
spp537<-function(a,b,c,d,e) {(c*e)-(e/a+d/a)^2+(0*(a+b+c+d+e))}
spp538<-function(a,b,c,d,e) {(c*e)-(e/b+e/a)^2+(0*(a+b+c+d+e))}
spp539<-function(a,b,c,d,e) {(c*e)-(2)^2+(0*(a+b+c+d+e))}
spp540<-function(a,b,c,d,e) {(c*e)-(20)^2+(0*(a+b+c+d+e))}

spp541<-function(a,b,c,d,e) {(c*d)-(e/a-b/a)^2 +(0*(a+b+c+d+e))}
spp542<-function(a,b,c,d,e) {(c*d)-(e-c)^2+(0*(a+b+c+d+e))}
spp543<-function(a,b,c,d,e) {(c*d)-(e-d)^2+(0*(a+b+c+d+e))}
spp544<-function(a,b,c,d,e) {(c*d)-(e-a)^2+(0*(a+b+c+d+e))}
spp545<-function(a,b,c,d,e) {(c*d)-(e/a+b/a)^2+(0*(a+b+c+d+e))}
spp546<-function(a,b,c,d,e) {(c*d)-(e/a+c/a)^2+(0*(a+b+c+d+e))}
spp547<-function(a,b,c,d,e) {(c*d)-(e/a+d/a)^2+(0*(a+b+c+d+e))}
spp548<-function(a,b,c,d,e) {(c*d)-(e/a+e/a)^2+(0*(a+b+c+d+e))}
spp549<-function(a,b,c,d,e) {(c*d)-(2)^2+(0*(a+b+c+d+e))}
spp550<-function(a,b,c,d,e) {(c*d)-(20)^2+(0*(a+b+c+d+e))}

spp551<-function(a,b,c,d,e) {(b*e)-(e-b)^2 +(0*(a+b+c+d+e))}
spp552<-function(a,b,c,d,e) {(b*e)-(e-c)^2+(0*(a+b+c+d+e))}
spp553<-function(a,b,c,d,e) {(b*e)-(e-d)^2+(0*(a+b+c+d+e))}
spp554<-function(a,b,c,d,e) {(b*e)-(e-a)^2+(0*(a+b+c+d+e))}
spp555<-function(a,b,c,d,e) {(b*e)-(e/c+b/c)^2+(0*(a+b+c+d+e))}
spp556<-function(a,b,c,d,e) {(b*e)-(e+c)^2+(0*(a+b+c+d+e))}
spp557<-function(a,b,c,d,e) {(b*e)-(e+d)^2+(0*(a+b+c+d+e))}
spp558<-function(a,b,c,d,e) {(b*e)-(e+e)^2+(0*(a+b+c+d+e))}
spp559<-function(a,b,c,d,e) {(b*e)-(2)^2+(0*(a+b+c+d+e))}
spp560<-function(a,b,c,d,e) {(b*e)-(20)^2+(0*(a+b+c+d+e))}

spp561<-function(a,b,c,d,e) {(a*e)-(e-b)^2 +(0*(a+b+c+d+e))}
spp562<-function(a,b,c,d,e) {(a*e)-(e-c)^2+(0*(a+b+c+d+e))}
spp563<-function(a,b,c,d,e) {(a*e)-(e-d)^2+(0*(a+b+c+d+e))}
spp564<-function(a,b,c,d,e) {(a*e)-(e-a)^2+(0*(a+b+c+d+e))}
spp565<-function(a,b,c,d,e) {(a*e)-(e/c+b/c)^2+(0*(a+b+c+d+e))}
spp566<-function(a,b,c,d,e) {(a*e)-(e+c)^2+(0*(a+b+c+d+e))}
spp567<-function(a,b,c,d,e) {(a*e)-(e+d)^2+(0*(a+b+c+d+e))}
spp568<-function(a,b,c,d,e) {(a*e)-(e+e)^2+(0*(a+b+c+d+e))}
spp569<-function(a,b,c,d,e) {(a*e)-(2)^2+(0*(a+b+c+d+e))}
spp570<-function(a,b,c,d,e) {(a*e)-(20)^2+(0*(a+b+c+d+e))}

spp571<-function(a,b,c,d,e) {(a*c)-(e/d-b/d)^2 +(0*(a+b+c+d+e))}
spp572<-function(a,b,c,d,e) {(a*c)-(e-c)^2+(0*(a+b+c+d+e))}
spp573<-function(a,b,c,d,e) {(a*c)-(e-d)^2+(0*(a+b+c+d+e))}
spp574<-function(a,b,c,d,e) {(a*c)-(e-a)^2+(0*(a+b+c+d+e))}
spp575<-function(a,b,c,d,e) {(a*c)-(e/d+b/d)^2+(0*(a+b+c+d+e))}
spp576<-function(a,b,c,d,e) {(a*c)-(e+c)^2+(0*(a+b+c+d+e))}
spp577<-function(a,b,c,d,e) {(a*c)-(e+d)^2+(0*(a+b+c+d+e))}
spp578<-function(a,b,c,d,e) {(a*c)-(e+e)^2+(0*(a+b+c+d+e))}
spp579<-function(a,b,c,d,e) {(a*c)-(2)^2+(0*(a+b+c+d+e))}
spp580<-function(a,b,c,d,e) {(a*c)-(20)^2+(0*(a+b+c+d+e))}

spp581<-function(a,b,c,d,e) {(b*c)-(e-b)^2 +(0*(a+b+c+d+e))}
spp582<-function(a,b,c,d,e) {(b*c)-(e-c)^2+(0*(a+b+c+d+e))}
spp583<-function(a,b,c,d,e) {(b*c)-(e-d)^2+(0*(a+b+c+d+e))}
spp584<-function(a,b,c,d,e) {(b*c)-(e-a)^2+(0*(a+b+c+d+e))}
spp585<-function(a,b,c,d,e) {(b*c)-(e/a+b/a)^2+(0*(a+b+c+d+e))}
spp586<-function(a,b,c,d,e) {(b*c)-(e+c)^2+(0*(a+b+c+d+e))}
spp587<-function(a,b,c,d,e) {(b*c)-(e+d)^2+(0*(a+b+c+d+e))}
spp588<-function(a,b,c,d,e) {(b*c)-(e+e)^2+(0*(a+b+c+d+e))}
spp589<-function(a,b,c,d,e) {(b*c)-(2)^2+(0*(a+b+c+d+e))}
spp590<-function(a,b,c,d,e) {(b*c)-(20)^2+(0*(a+b+c+d+e))}

spp591<-function(a,b,c,d,e) {(d*a)-(e/c-b/c)^2 +(0*(a+b+c+d+e))}
spp592<-function(a,b,c,d,e) {(d*a)-(e-c)^2+(0*(a+b+c+d+e))}
spp593<-function(a,b,c,d,e) {(d*a)-(e-d)^2+(0*(a+b+c+d+e))}
spp594<-function(a,b,c,d,e) {(d*a)-(e-a)^2+(0*(a+b+c+d+e))}
spp595<-function(a,b,c,d,e) {(d*a)-(e/c+b/c)^2+(0*(a+b+c+d+e))}
spp596<-function(a,b,c,d,e) {(d*a)-(e+c)^2+(0*(a+b+c+d+e))}
spp597<-function(a,b,c,d,e) {(d*a)-(e+d)^2+(0*(a+b+c+d+e))}
spp598<-function(a,b,c,d,e) {(d*a)-(e+e)^2+(0*(a+b+c+d+e))}
spp599<-function(a,b,c,d,e) {(d*a)-(2)^2+(0*(a+b+c+d+e))}
spp600<-function(a,b,c,d,e) {(d*a)-(20)^2+(0*(a+b+c+d+e))}######################################

spp601<-function(a,b,c,d,e) {(c*d*e)-(a)^2 +(0*(a+b+c+d+e))}
spp602<-function(a,b,c,d,e) {(c*d*e)-(b)^2+(0*(a+b+c+d+e))}
spp603<-function(a,b,c,d,e) {(c*d*e)-(c)^2+(0*(a+b+c+d+e))}
spp604<-function(a,b,c,d,e) {(c*d*e)-(d)^2+(0*(a+b+c+d+e))}
spp605<-function(a,b,c,d,e) {(c*d*e)-(e)^2+(0*(a+b+c+d+e))}
spp606<-function(a,b,c,d,e) {(c*d*e)-(2*a)^2+(0*(a+b+c+d+e))}
spp607<-function(a,b,c,d,e) {(c*d*e)-(2*b)^2+(0*(a+b+c+d+e))}
spp608<-function(a,b,c,d,e) {(c*d*e)-(2*c)^2+(0*(a+b+c+d+e))}
spp609<-function(a,b,c,d,e) {(c*d*e)-(2*d)^2+(0*(a+b+c+d+e))}
spp610<-function(a,b,c,d,e) {(c*d*e)-(2*e)^2+(0*(a+b+c+d+e))}

spp611<-function(a,b,c,d,e) {(d*e)-(a)^2 +(0*(a+b+c+d+e))}
spp612<-function(a,b,c,d,e) {(d*e)-(b)^2+(0*(a+b+c+d+e))}
spp613<-function(a,b,c,d,e) {(d*e)-(c)^2+(0*(a+b+c+d+e))}
spp614<-function(a,b,c,d,e) {(d*e)-(d)^2+(0*(a+b+c+d+e))}
spp615<-function(a,b,c,d,e) {(d*e)-(e)^2+(0*(a+b+c+d+e))}
spp616<-function(a,b,c,d,e) {(d*e)-(2*a)^2+(0*(a+b+c+d+e))}
spp617<-function(a,b,c,d,e) {(d*e)-(2*b)^2+(0*(a+b+c+d+e))}
spp618<-function(a,b,c,d,e) {(d*e)-(2*c)^2+(0*(a+b+c+d+e))}
spp619<-function(a,b,c,d,e) {(d*e)-(2*d)^2+(0*(a+b+c+d+e))}
spp620<-function(a,b,c,d,e) {(d*e)-(2*e/a)^2+(0*(a+b+c+d+e))}

spp621<-function(a,b,c,d,e) {(a*b)-(a)^2 +(0*(a+b+c+d+e))}
spp622<-function(a,b,c,d,e) {(a*b)-(b)^2+(0*(a+b+c+d+e))}
spp623<-function(a,b,c,d,e) {(a*b)-(c)^2+(0*(a+b+c+d+e))}
spp624<-function(a,b,c,d,e) {(a*b)-(d)^2+(0*(a+b+c+d+e))}
spp625<-function(a,b,c,d,e) {(a*b)-(e)^2+(0*(a+b+c+d+e))}
spp626<-function(a,b,c,d,e) {(a*b)-(2*a)^2+(0*(a+b+c+d+e))}
spp627<-function(a,b,c,d,e) {(a*b)-(2*b/c)^2+(0*(a+b+c+d+e))}
spp628<-function(a,b,c,d,e) {(a*b)-(2*c)^2+(0*(a+b+c+d+e))}
spp629<-function(a,b,c,d,e) {(a*b)-(2*d)^2+(0*(a+b+c+d+e))}
spp630<-function(a,b,c,d,e) {(a*b)-(2*e)^2+(0*(a+b+c+d+e))}

spp631<-function(a,b,c,d,e) {(c*e)-(a)^2 +(0*(a+b+c+d+e))}
spp632<-function(a,b,c,d,e) {(c*e)-(b)^2+(0*(a+b+c+d+e))}
spp633<-function(a,b,c,d,e) {(c*e)-(c)^2+(0*(a+b+c+d+e))}
spp634<-function(a,b,c,d,e) {(c*e)-(d)^2+(0*(a+b+c+d+e))}
spp635<-function(a,b,c,d,e) {(c*e)-(e)^2+(0*(a+b+c+d+e))}
spp636<-function(a,b,c,d,e) {(c*e)-(2*a)^2+(0*(a+b+c+d+e))}
spp637<-function(a,b,c,d,e) {(c*e)-(2*b)^2+(0*(a+b+c+d+e))}
spp638<-function(a,b,c,d,e) {(c*e)-(2*c)^2+(0*(a+b+c+d+e))}
spp639<-function(a,b,c,d,e) {(c*e)-(2*d)^2+(0*(a+b+c+d+e))}
spp640<-function(a,b,c,d,e) {(c*e)-(2*e/d)^2+(0*(a+b+c+d+e))}

spp641<-function(a,b,c,d,e) {(c*d)-(a)^2 +(0*(a+b+c+d+e))}
spp642<-function(a,b,c,d,e) {(c*d)-(b)^2+(0*(a+b+c+d+e))}
spp643<-function(a,b,c,d,e) {(c*d)-(c)^2+(0*(a+b+c+d+e))}
spp644<-function(a,b,c,d,e) {(c*d)-(d)^2+(0*(a+b+c+d+e))}
spp645<-function(a,b,c,d,e) {(c*d)-(e)^2+(0*(a+b+c+d+e))}
spp646<-function(a,b,c,d,e) {(c*d)-(2*a)^2+(0*(a+b+c+d+e))}
spp647<-function(a,b,c,d,e) {(c*d)-(3*b)^2+(0*(a+b+c+d+e))}
spp648<-function(a,b,c,d,e) {(c*d)-(2*c/a)^2+(0*(a+b+c+d+e))}
spp649<-function(a,b,c,d,e) {(c*d)-(2*d/a)^2+(0*(a+b+c+d+e))}
spp650<-function(a,b,c,d,e) {(c*d)-(2*e/a)^2+(0*(a+b+c+d+e))}

spp651<-function(a,b,c,d,e) {(b*e)-(a)^2 +(0*(a+b+c+d+e))}
spp652<-function(a,b,c,d,e) {(b*e)-(b)^2+(0*(a+b+c+d+e))}
spp653<-function(a,b,c,d,e) {(b*e)-(c)^2+(0*(a+b+c+d+e))}
spp654<-function(a,b,c,d,e) {(b*e)-(d)^2+(0*(a+b+c+d+e))}
spp655<-function(a,b,c,d,e) {(b*e)-(e)^2+(0*(a+b+c+d+e))}
spp656<-function(a,b,c,d,e) {(b*e)-(2*a)^2+(0*(a+b+c+d+e))}
spp657<-function(a,b,c,d,e) {(b*e)-(2*b)^2+(0*(a+b+c+d+e))}
spp658<-function(a,b,c,d,e) {(b*e)-(2*c)^2+(0*(a+b+c+d+e))}
spp659<-function(a,b,c,d,e) {(b*e)-(2*d)^2+(0*(a+b+c+d+e))}
spp660<-function(a,b,c,d,e) {(b*e)-(2*e)^2+(0*(a+b+c+d+e))}

spp661<-function(a,b,c,d,e) {(a*e)-(a)^2 +(0*(a+b+c+d+e))}
spp662<-function(a,b,c,d,e) {(a*e)-(b)^2+(0*(a+b+c+d+e))}
spp663<-function(a,b,c,d,e) {(a*e)-(c)^2+(0*(a+b+c+d+e))}
spp664<-function(a,b,c,d,e) {(a*e)-(d)^2+(0*(a+b+c+d+e))}
spp665<-function(a,b,c,d,e) {(a*e)-(e)^2+(0*(a+b+c+d+e))}
spp666<-function(a,b,c,d,e) {(a*e)-(2*a)^2+(0*(a+b+c+d+e))}
spp667<-function(a,b,c,d,e) {(a*e)-(2*b)^2+(0*(a+b+c+d+e))}
spp668<-function(a,b,c,d,e) {(a*e)-(2*c)^2+(0*(a+b+c+d+e))}
spp669<-function(a,b,c,d,e) {(a*e)-(2*d)^2+(0*(a+b+c+d+e))}
spp670<-function(a,b,c,d,e) {(a*e)-(2*e)^2+(0*(a+b+c+d+e))}

spp671<-function(a,b,c,d,e) {(a*c)-(a)^2 +(0*(a+b+c+d+e))}
spp672<-function(a,b,c,d,e) {(a*c)-(b)^2+(0*(a+b+c+d+e))}
spp673<-function(a,b,c,d,e) {(a*c)-(c)^2+(0*(a+b+c+d+e))}
spp674<-function(a,b,c,d,e) {(a*c)-(d)^2+(0*(a+b+c+d+e))}
spp675<-function(a,b,c,d,e) {(a*c)-(e)^2+(0*(a+b+c+d+e))}
spp676<-function(a,b,c,d,e) {(a*c)-(2*a/d)^2+(0*(a+b+c+d+e))}
spp677<-function(a,b,c,d,e) {(a*c)-(2*b)^2+(0*(a+b+c+d+e))}
spp678<-function(a,b,c,d,e) {(a*c)-(2*c)^2+(0*(a+b+c+d+e))}
spp679<-function(a,b,c,d,e) {(a*c)-(2*d)^2+(0*(a+b+c+d+e))}
spp680<-function(a,b,c,d,e) {(a*c)-(2*e)^2+(0*(a+b+c+d+e))}

spp681<-function(a,b,c,d,e) {(b*c)-(a)^2 +(0*(a+b+c+d+e))}
spp682<-function(a,b,c,d,e) {(b*c)-(b)^2+(0*(a+b+c+d+e))}
spp683<-function(a,b,c,d,e) {(b*c)-(c)^2+(0*(a+b+c+d+e))}
spp684<-function(a,b,c,d,e) {(b*c)-(d)^2+(0*(a+b+c+d+e))}
spp685<-function(a,b,c,d,e) {(b*c)-(e)^2+(0*(a+b+c+d+e))}
spp686<-function(a,b,c,d,e) {(b*c)-(2*a)^2+(0*(a+b+c+d+e))}
spp687<-function(a,b,c,d,e) {(b*c)-(2*b/c)^2+(0*(a+b+c+d+e))}
spp688<-function(a,b,c,d,e) {(b*c)-(2*c)^2+(0*(a+b+c+d+e))}
spp689<-function(a,b,c,d,e) {(b*c)-(2*d)^2+(0*(a+b+c+d+e))}
spp690<-function(a,b,c,d,e) {(b*c)-(2*e)^2+(0*(a+b+c+d+e))}

spp691<-function(a,b,c,d,e) {(d*a)-(a)^2 +(0*(a+b+c+d+e))}
spp692<-function(a,b,c,d,e) {(d*a)-(b)^2+(0*(a+b+c+d+e))}
spp693<-function(a,b,c,d,e) {(d*a)-(c)^2+(0*(a+b+c+d+e))}
spp694<-function(a,b,c,d,e) {(d*a)-(d)^2+(0*(a+b+c+d+e))}
spp695<-function(a,b,c,d,e) {(d*a)-(e)^2+(0*(a+b+c+d+e))}
spp696<-function(a,b,c,d,e) {(d*a)-(2*a/c)^2+(0*(a+b+c+d+e))}
spp697<-function(a,b,c,d,e) {(d*a)-(2*b/c)^2+(0*(a+b+c+d+e))}
spp698<-function(a,b,c,d,e) {(d*a)-(2*c)^2+(0*(a+b+c+d+e))}
spp699<-function(a,b,c,d,e) {(d*a)-(2*d)^2+(0*(a+b+c+d+e))}
spp700<-function(a,b,c,d,e) {(d*a)-(2*e)^2+(0*(a+b+c+d+e))}

################################




# wrapper for making spp list: ####
make.spList<-function(n, replace=F){
  AllSpp<-c(paste0("spp", c(1:1724), sep=""))
  AllSpp<-lapply(AllSpp, get)
  AllSpp<-unlist(AllSpp)
  
  C<-base::sample(AllSpp, n, replace)
  C
}
# all

# generate the community workflow ####
AllSpp<-c(paste0("spp", c(1:1724), sep="")) # trying to make a quick list of all functions
AllSpp<-lapply(AllSpp, get)
AllSpp<-unlist(AllSpp)

seeds<-round(rnorm(50,5000,100)) 
set.seed(seeds[1]) # 4945
Comm1<-base::sample(AllSpp, 200, replace=F)
set.seed(seeds[2]) # 4902
Comm2<-base::sample(AllSpp, 200, replace=F)
set.seed(seeds[4]) # 5098
Comm3<-base::sample(AllSpp, 200, replace=F)

names(Comm1)
names(Comm2)
names(Comm3)




library(reshape2)
f1c1<-c(5,5,5,5,5,5)
f1c2<-c(1,3,10,15,3,15)
f1c3<-c(0.5,0.5,3,3,1,5)
F1.frame<-mapply(rnorm, f1c1,f1c2,f1c3)
F1<-melt(F1.frame)

#F2
f2c1<-c(5,5,5,5,5,5)
f2c2<-c(34,30,50,55,35,60)
f2c3<-c(0.5,0.5,3,3,1,5)
F2.frame<-mapply(rnorm, f2c1,f2c2,f2c3)
F2<-melt(F2.frame)

#F3
f3c1<-c(5,5,5,5,5,5)
f3c2<-c(1,3,10,15,3,15)
f3c3<-c(0.5,0.5,3,3,1,5)
F3.frame<-mapply(rnorm, f3c1,f3c2,f3c3)
F3<-melt(F3.frame)

#F4
f4c1<-c(5,5,5,5,5,5)
f4c2<-c(1,3,10,15,3,15)
f4c3<-c(0.5,0.5,3,3,1,5)
F4.frame<-mapply(rnorm, f4c1,f4c2,f4c3)
F4<-melt(F4.frame)

#F5
f5c1<-c(5,5,5,5,5,5)
f5c2<-c(1,3,10,15,3,15)
f5c3<-c(0.5,0.5,3,3,1,5)
F5.frame<-mapply(rnorm, f5c1,f5c2,f5c3)
F5<-melt(F5.frame)

Factors<-data.frame(F1$value,F2$value,F3$value,F4$value,F5$value) # environment
Sites<-c(paste0("Site", 1:30))
rownames(Factors)<-Sites
colnames(Factors)<-c("F1","F2","F3","F4","F5")
head(Factors)

saveRDS(Factors, "~/Documents/Github/ModelMicrobiome/Model_environment.RDS") # save environmental gradient!!

#### output response table ###

make.comm<-function(Comm1, Factors){
otu<-matrix(data=NA, nrow=nrow(Factors), ncol = length(Comm1))
Sites<-c(paste0("Site", 1:30))
for(i in 1:length(Comm1)) {
  for(row in 1:nrow(Factors)){
   otu[row,i]<-do.call(Comm1[[i]], list(Factors[row,1],Factors[row,2],Factors[row,3],Factors[row,4],Factors[row,5]))
      }
}

#otu<-make.comm(Comm1, Factors)
row.names(otu)<-Sites
colnames(otu)<-names(Comm1)
otu[otu<0]<-0
otu<-round(otu)
otu<-otu_table(otu, taxa_are_rows = FALSE)
Sa<-sample_data(Factors)
out<-phyloseq(otu, Sa)
out}


t.ab<-sample_sums(out)
sample_data(out)$total_abund<-t.ab
saveRDS(out, "~/Documents/GitHub/ModelMicrobiome/PS_envComm.RDS")
#####
{
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function is meaningful only for integers (counts)")
  x <- as.matrix(x)
  if (ncol(x) == 1) 
    x <- t(x)
  if (length(sample) > 1 && length(sample) != nrow(x)) 
    stop(gettextf("length of 'sample' and number of rows of 'x' do not match"))
  sample <- rep(sample, length = nrow(x))
  colnames(x) <- colnames(x, do.NULL = FALSE)
  nm <- colnames(x)
  if (any(rowSums(x) < sample)) 
    warning("Some row sums < 'sample' and are not rarefied")
  for (i in 1:nrow(x)) {
    if (sum(x[i, ]) <= sample[i]) 
      next 
    row <- sample(rep(nm, times = x[i, ]), sample[i], replace = T)
    row <- table(row)
    ind <- names(row)
    x[i, ] <- 0
    x[i, ind] <- row
    }
  }
  x
}

rrarefy2<-function (x, sample, replace, prob=NULL) 
{
    x <- as.matrix(x)
    if (!identical(all.equal(x, round(x)), TRUE)) 
        stop("function is meaningful only for integers (counts)")
    if (!is.integer(x)) 
        x <- round(x)
    if (ncol(x) == 1) 
        x <- t(x)
    if (length(sample) > 1 && length(sample) != nrow(x)) 
        stop(gettextf("length of 'sample' and number of rows of 'x' do not match"))
    sample <- rep(sample, length = nrow(x))
    if (any(rowSums(x) < sample)) 
        warning("some row sums < 'sample' and are not rarefied")
    for (i in 1:nrow(x)) {
        x[i, ] <- .Call(do_rrarefy, x[i, ], sample[i])
    }
    x
}
#####


#### Select core 10 taxa found at all sites
# use lists
group<-c(1:5)

df<-data.frame("group"=numeric(), "spp"=numeric())
for (g in 2:group) {
a<-names(sample(AllSpp, 10, replace=F))
b<-c(rep(g, 10))

df <- rbind(df, data.frame(b,a))

}

# need data frame with structure (sites + groups)
library(plyr)

global.spp<-names(sample(AllSpp, 10, replace=F))

group.spp<-NULL
group.spp$group1<-names(sample(AllSpp, 20, replace=F))
group.spp$group2<-names(sample(AllSpp, 20, replace=F))
group.spp$group3<-names(sample(AllSpp, 20, replace=F))
group.spp$group4<-names(sample(AllSpp, 20, replace=F))
group.spp$group5<-names(sample(AllSpp, 20, replace=F))
group.spp$group6<-names(sample(AllSpp, 20, replace=F))

rando.spp<-NULL
rando.spp$Site1<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group1), global.spp))
rando.spp$Site2<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group1), global.spp))
rando.spp$Site3<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group1), global.spp))
rando.spp$Site4<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group1), global.spp))
rando.spp$Site5<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group1), global.spp))
rando.spp$Site6<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group2), global.spp))
rando.spp$Site7<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group2), global.spp))
rando.spp$Site8<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group2), global.spp))
rando.spp$Site9<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group2), global.spp))
rando.spp$Site10<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group2), global.spp))
rando.spp$Site11<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group3), global.spp))
rando.spp$Site12<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group3), global.spp))
rando.spp$Site13<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group3), global.spp))
rando.spp$Site14<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group3), global.spp))
rando.spp$Site15<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group3), global.spp))
rando.spp$Site16<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group4), global.spp))
rando.spp$Site17<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group4), global.spp))
rando.spp$Site18<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group4), global.spp))
rando.spp$Site19<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group4), global.spp))
rando.spp$Site20<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group4), global.spp))
rando.spp$Site21<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group5), global.spp))
rando.spp$Site22<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group5), global.spp))
rando.spp$Site23<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group5), global.spp))
rando.spp$Site24<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group5), global.spp))
rando.spp$Site25<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group5), global.spp))
rando.spp$Site26<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group6), global.spp))
rando.spp$Site27<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group6), global.spp))
rando.spp$Site28<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group6), global.spp))
rando.spp$Site29<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group6), global.spp))
rando.spp$Site30<-unique(c(names(sample(AllSpp, 70, replace=F)), c(group.spp$group6), global.spp))

rando.spp
final.list<-NULL
final.list$Site1<-lapply(rando.spp$Site1, get)
final.list$Site2<-lapply(rando.spp$Site2, get)
final.list$Site3<-lapply(rando.spp$Site3, get)
final.list$Site4<-lapply(rando.spp$Site4, get)
final.list$Site5<-lapply(rando.spp$Site5, get)
final.list$Site6<-lapply(rando.spp$Site6, get)
final.list$Site7<-lapply(rando.spp$Site7, get)
final.list$Site8<-lapply(rando.spp$Site8, get)
final.list$Site9<-lapply(rando.spp$Site9, get)
final.list$Site10<-lapply(rando.spp$Site10, get)
final.list$Site11<-lapply(rando.spp$Site11, get)
final.list$Site12<-lapply(rando.spp$Site12, get)
final.list$Site13<-lapply(rando.spp$Site13, get)
final.list$Site14<-lapply(rando.spp$Site14, get)
final.list$Site15<-lapply(rando.spp$Site15, get)
final.list$Site16<-lapply(rando.spp$Site16, get)
final.list$Site17<-lapply(rando.spp$Site17, get)
final.list$Site18<-lapply(rando.spp$Site18, get)
final.list$Site19<-lapply(rando.spp$Site19, get)
final.list$Site20<-lapply(rando.spp$Site20, get)
final.list$Site21<-lapply(rando.spp$Site21, get)
final.list$Site22<-lapply(rando.spp$Site22, get)
final.list$Site23<-lapply(rando.spp$Site23, get)
final.list$Site24<-lapply(rando.spp$Site24, get)
final.list$Site25<-lapply(rando.spp$Site25, get)
final.list$Site26<-lapply(rando.spp$Site26, get)
final.list$Site27<-lapply(rando.spp$Site27, get)
final.list$Site28<-lapply(rando.spp$Site28, get)
final.list$Site29<-lapply(rando.spp$Site29, get)
final.list$Site30<-lapply(rando.spp$Site30, get)

names(final.list)<-names(rando.spp)



make.comm2<-function(rando.spp, Factors){
l1<-NULL
for (i in 1:length(rando.spp[[1]])){
  l1[i]<-do.call(rando.spp[[1]][i], list(Factors[1,1],Factors[1,2],Factors[1,3],Factors[1,4],Factors[1,5]))
  }
#l1<-data.frame("Site1"=l1, "Spp"=rando.spp[[1]])
names(l1)<-rando.spp[[1]]
for (r in 2:nrow(Factors)) # for each site...
{ l2<-NULL
  for (i in 1:length(rando.spp[[r]])){  # for each species in site...
    l2[i]<-do.call(rando.spp[[r]][i], list(Factors[r,1],Factors[r,2],Factors[r,3],Factors[r,4],Factors[r,5]))
    }
  names(l2)<-rando.spp[[r]]
  l1<-merge(as.data.frame(l1),as.data.frame(l2), by=0, all=T)
  rownames(l1)<-l1$Row.names
  colnames(l1)[colnames(l1) == "l1"] <- "Site1"
 colnames(l1)[colnames(l1) == "l2"] <- paste("Site", r, sep="")
 l1<-l1[,-1]
  }
l1<-round(l1)
l1[is.na(l1)]<-0
l1[l1<0]<-0
l1
}

p.css<-function(s, ...){
  require(phyloseq)
  t<-subset_samples(..., Factor=="one"|Factor==s)
  t
}

Make.Css<-function(ps, list){
  require(phyloseq)
  require(vegan)
  require(metagenomeSeq)
  
  l1<-lapply(list, p.css)
  l2<-ldply(l1, function(ps){
  a<-phyloseq_to_metagenomeSeq(ps)
  b<-cumNorm(a, p=0.75)
  sam<-pData(b)
  mod <-  model.matrix(~1+Factor, data = sam)
  res<-fitFeatureModel(b, mod)
  res@pvalues})
  
  logFC <- MRcoefs(res, number = nrow(b))
  feature_names <- rownames(logFC[which(logFC[which(abs(logFC$logFC) > 1),]$adjPvalues < .1),]) 
  fSelection <- generateSelection(feature_names = feature_names , aggregation_level = aggregation_level, selection_type =2)
  MRcounts(b, norm = T)
  
}


anova.test<-function(ps,dependent){
  require(phyloseq)
  out<-matrix(nrow=ntaxa(ps), ncol=2)
  for(i in 1:ntaxa(ps)){
    a<-aov()
    out[i, 1]<-a$
  }
}