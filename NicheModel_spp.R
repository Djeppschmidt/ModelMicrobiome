#define species functions
#species response (substrate/growth function)+(environmental responces)
#define boundary ranges; the functions must produce at least some positive numbers across the range of input values

# initial model
# Linear models

# high abundance (1000s - beyond)
spp1<-function(a,b,c,d,e) {(0.001(a-50)^3+3)+(10*b)+(-0.1(c-50)^2+50)+d+e}
spp2<-function(a,b,c,d,e) {(0.001(a-50)^3+3)*(10*b)+(-0.1(c-50)^2+50)+d+e}
spp3<-function(a,b,c,d,e) {((1/e)*(a-50)^3+3)+(10*b)+(-0.1(c-50)^2+50)+d}
spp4<-function(a,b,c,d,e) {((1/e)*(a-50)^3+3)+(10*b)+(-(1/d)*(c-50)^2+50)}
spp5<-function(a,b,c,d,e) {((1/e)^2*(a-50)^3+3)+(10*b)+(-(1/d)^2*(c-50)^2+50)}
spp6<-function(a,b,c,d,e) {}
spp7<-function(a,b,c,d,e) {}
spp8<-function(a,b,c,d,e) {}
spp9<-function(a,b,c,d,e) {}
spp10<-function(a,b,c,d,e) {}
spp11<-function(a,b,c,d,e) {}
spp12<-function(a,b,c,d,e) {}
spp13<-function(a,b,c,d,e) {}
spp14<-function(a,b,c,d,e) {}
spp15<-function(a,b,c,d,e) {}
spp16<-function(a,b,c,d,e) {}
spp17<-function(a,b,c,d,e) {}
spp18<-function(a,b,c,d,e) {}
spp19<-function(a,b,c,d,e) {}
spp20<-function(a,b,c,d,e) {}
spp21<-function(a,b,c,d,e) {}
spp22<-function(a,b,c,d,e) {}
spp23<-function(a,b,c,d,e) {}
spp24<-function(a,b,c,d,e) {}
spp25<-function(a,b,c,d,e) {}
spp26<-function(a,b,c,d,e) {}
spp27<-function(a,b,c,d,e) {}
spp28<-function(a,b,c,d,e) {}
spp29<-function(a,b,c,d,e) {}
spp30<-function(a,b,c,d,e) {}
spp31<-function(a,b,c,d,e) {}

#medium abundance (100s - 1000s)

spp32<-function(a,b,c,d,e) {}
spp33<-function(a,b,c,d,e) {}
spp34<-function(a,b,c,d,e) {}
spp35<-function(a,b,c,d,e) {}
spp36<-function(a,b,c,d,e) {}
spp37<-function(a,b,c,d,e) {}
spp38<-function(a,b,c,d,e) {}
spp39<-function(a,b,c,d,e) {}
spp40<-function(a,b,c,d,e) {}
spp41<-function(a,b,c,d,e) {}
spp42<-function(a,b,c,d,e) {}
spp43<-function(a,b,c,d,e) {}
spp44<-function(a,b,c,d,e) {}
spp45<-function(a,b,c,d,e) {}
spp46<-function(a,b,c,d,e) {}
spp47<-function(a,b,c,d,e) {}
spp48<-function(a,b,c,d,e) {}
spp49<-function(a,b,c,d,e) {}
spp50<-function(a,b,c,d,e) {}
spp51<-function(a,b,c,d,e) {}
spp52<-function(a,b,c,d,e) {}
spp53<-function(a,b,c,d,e) {}
spp54<-function(a,b,c,d,e) {}
spp55<-function(a,b,c,d,e) {}
spp56<-function(a,b,c,d,e) {}
spp57<-function(a,b,c,d,e) {}
spp58<-function(a,b,c,d,e) {}
spp59<-function(a,b,c,d,e) {}
spp60<-function(a,b,c,d,e) {}
spp61<-function(a,b,c,d,e) {}
spp62<-function(a,b,c,d,e) {}

#low abundances (1-100)

spp63<-function(a,b,c,d,e) {}
spp64<-function(a,b,c,d,e) {}
spp65<-function(a,b,c,d,e) {}
spp66<-function(a,b,c,d,e) {}
spp67<-function(a,b,c,d,e) {}
spp68<-function(a,b,c,d,e) {}
spp69<-function(a,b,c,d,e) {}
spp70<-function(a,b,c,d,e) {}
spp71<-function(a,b,c,d,e) {}
spp72<-function(a,b,c,d,e) {}
spp73<-function(a,b,c,d,e) {}
spp74<-function(a,b,c,d,e) {}
spp75<-function(a,b,c,d,e) {}
spp76<-function(a,b,c,d,e) {}
spp77<-function(a,b,c,d,e) {}
spp78<-function(a,b,c,d,e) {}
spp79<-function(a,b,c,d,e) {}
spp80<-function(a,b,c,d,e) {}
spp81<-function(a,b,c,d,e) {}
spp82<-function(a,b,c,d,e) {}
spp83<-function(a,b,c,d,e) {}
spp84<-function(a,b,c,d,e) {}
spp85<-function(a,b,c,d,e) {}
spp86<-function(a,b,c,d,e) {}
spp87<-function(a,b,c,d,e) {}
spp88<-function(a,b,c,d,e) {}
spp89<-function(a,b,c,d,e) {}
spp90<-function(a,b,c,d,e) {}
spp91<-function(a,b,c,d,e) {}
spp92<-function(a,b,c,d,e) {}
spp93<-function(a,b,c,d,e) {}
spp94<-function(a,b,c,d,e) {}
spp95<-function(a,b,c,d,e) {}
spp96<-function(a,b,c,d,e) {}
spp97<-function(a,b,c,d,e) {}
spp98<-function(a,b,c,d,e) {}
spp99<-function(a,b,c,d,e) {}
spp100<-function(a,b,c,d,e) {}


####
####
#### Non-linear models
####
####

## Low
# multinomial/exp
spp101<-function(a,b,c,d,e) {(0.001(a-50)^3+3)+(10*b)+(-0.1(c-50)^2+50)+d+e}
spp102<-function(a,b,c,d,e) {(0.001(a-50)^3+3)*(10*b)+(-0.1(c-50)^2+50)+d+e}
spp103<-function(a,b,c,d,e) {((1/e)*(a-50)^3+3)+(10*b)+(-0.1(c-50)^2+50)+d}
spp104<-function(a,b,c,d,e) {((1/e)*(a-50)^3+3)+(10*b)+(-(1/d)*(c-50)^2+50)}
spp105<-function(a,b,c,d,e) {((1/e)^2*(a-50)^3+3)+(10*b)+(-(1/d)^2*(c-50)^2+50)}
spp106<-function(a,b,c,d,e) {((1/d)*(a-50)^3+3)+(10*b)+(-(1/e)*(c-50)^2+50)}
spp107<-function(a,b,c,d,e) {((1/e)*(a-100)^3+3)*(10*b)+(-(1/d)*(c-50)^2+50)}
spp108<-function(a,b,c,d,e) {((1/e)^2*(a-100)^3+3)*(10*b)+(-(1/d)^2*(c-100)^2+50)}
spp109<-function(a,b,c,d,e) {(0.001(a-50)^3+3)*(10*b)*(1/d)+(-0.1(c-50)^2+50)+d+e}
spp110<-function(a,b,c,d,e) {(0.001(a-50)^3+3)*(10*b)*(-0.001(e-50)^2+50)+(-0.00000001(e-50)^5+10)+(-0.1(c-50)^2+50)+d}
# geometric
spp111<-function(a,b,c,d,e) {}
spp112<-function(a,b,c,d,e) {}
spp113<-function(a,b,c,d,e) {}
spp114<-function(a,b,c,d,e) {}
spp115<-function(a,b,c,d,e) {}
spp116<-function(a,b,c,d,e) {}
spp117<-function(a,b,c,d,e) {}
spp118<-function(a,b,c,d,e) {}
spp119<-function(a,b,c,d,e) {}
spp120<-function(a,b,c,d,e) {}
spp121<-function(a,b,c,d,e) {}

#log

spp122<-function(a,b,c,d,e) {}
spp123<-function(a,b,c,d,e) {}
spp124<-function(a,b,c,d,e) {}
spp125<-function(a,b,c,d,e) {}
spp126<-function(a,b,c,d,e) {}
spp127<-function(a,b,c,d,e) {}
spp128<-function(a,b,c,d,e) {}
spp129<-function(a,b,c,d,e) {}
spp130<-function(a,b,c,d,e) {}
spp131<-function(a,b,c,d,e) {}
spp132<-function(a,b,c,d,e) {}
spp133<-function(a,b,c,d,e) {}

# medium

spp134<-function(a,b,c,d,e) {}
spp135<-function(a,b,c,d,e) {}
spp136<-function(a,b,c,d,e) {}
spp137<-function(a,b,c,d,e) {}
spp138<-function(a,b,c,d,e) {}
spp139<-function(a,b,c,d,e) {}
spp140<-function(a,b,c,d,e) {}
spp141<-function(a,b,c,d,e) {}
spp142<-function(a,b,c,d,e) {}
spp143<-function(a,b,c,d,e) {}
spp144<-function(a,b,c,d,e) {}
spp145<-function(a,b,c,d,e) {}
spp146<-function(a,b,c,d,e) {}
spp147<-function(a,b,c,d,e) {}
spp148<-function(a,b,c,d,e) {}
spp149<-function(a,b,c,d,e) {}
spp150<-function(a,b,c,d,e) {}
spp151<-function(a,b,c,d,e) {}
spp152<-function(a,b,c,d,e) {}
spp153<-function(a,b,c,d,e) {}
spp154<-function(a,b,c,d,e) {}
spp155<-function(a,b,c,d,e) {}
spp156<-function(a,b,c,d,e) {}
spp157<-function(a,b,c,d,e) {}
spp158<-function(a,b,c,d,e) {}
spp159<-function(a,b,c,d,e) {}
spp160<-function(a,b,c,d,e) {}
spp161<-function(a,b,c,d,e) {}
spp162<-function(a,b,c,d,e) {}
spp163<-function(a,b,c,d,e) {}
spp164<-function(a,b,c,d,e) {}
spp165<-function(a,b,c,d,e) {}
spp166<-function(a,b,c,d,e) {}

# high

spp167<-function(a,b,c,d,e) {}
spp168<-function(a,b,c,d,e) {}
spp169<-function(a,b,c,d,e) {}
spp170<-function(a,b,c,d,e) {}
spp171<-function(a,b,c,d,e) {}
spp172<-function(a,b,c,d,e) {}
spp173<-function(a,b,c,d,e) {}
spp174<-function(a,b,c,d,e) {}
spp175<-function(a,b,c,d,e) {}
spp176<-function(a,b,c,d,e) {}
spp177<-function(a,b,c,d,e) {}
spp178<-function(a,b,c,d,e) {}
spp179<-function(a,b,c,d,e) {}
spp180<-function(a,b,c,d,e) {}
spp181<-function(a,b,c,d,e) {}
spp182<-function(a,b,c,d,e) {}
spp183<-function(a,b,c,d,e) {}
spp184<-function(a,b,c,d,e) {}
spp185<-function(a,b,c,d,e) {}
spp186<-function(a,b,c,d,e) {}
spp187<-function(a,b,c,d,e) {}
spp188<-function(a,b,c,d,e) {}
spp189<-function(a,b,c,d,e) {}
spp190<-function(a,b,c,d,e) {}
spp191<-function(a,b,c,d,e) {}
spp192<-function(a,b,c,d,e) {}
spp193<-function(a,b,c,d,e) {}
spp194<-function(a,b,c,d,e) {}
spp195<-function(a,b,c,d,e) {}
spp196<-function(a,b,c,d,e) {}
spp197<-function(a,b,c,d,e) {}
spp198<-function(a,b,c,d,e) {}
spp199<-function(a,b,c,d,e) {}
spp200<-function(a,b,c,d,e) {}



####
####
#### Random effects
####
####

# High
spp201<-function(a,b,c,d,e) {rnorm(1, 30000, 4000)+(0*(a+b+c+d+e))}
spp202<-function(a,b,c,d,e) {rnorm(1, 40000, 4000)+(0*(a+b+c+d+e))}
spp203<-function(a,b,c,d,e) {rnorm(1, 4000, 400)+(0*(a+b+c+d+e))}
spp204<-function(a,b,c,d,e) {rnorm(1, 4000, 40)+(0*(a+b+c+d+e))}
spp205<-function(a,b,c,d,e) {rnorm(1, 4000, 4000)+(0*(a+b+c+d+e))}
spp206<-function(a,b,c,d,e) {rnorm(1, 60000, 40000)+(0*(a+b+c+d+e))}
spp207<-function(a,b,c,d,e) {rnorm(1, 100000, 40000)+(0*(a+b+c+d+e))}
spp208<-function(a,b,c,d,e) {rnorm(1, 100000, 4000)+(0*(a+b+c+d+e))}
spp209<-function(a,b,c,d,e) {rnorm(1, 100000, 400)+(0*(a+b+c+d+e))}
spp210<-function(a,b,c,d,e) {rnorm(1, 100000, 40)+(0*(a+b+c+d+e))}
spp211<-function(a,b,c,d,e) {rnorm(1, 100000, 4)+(0*(a+b+c+d+e))}
spp212<-function(a,b,c,d,e) {rnorm(1, 10000, 40000)+(0*(a+b+c+d+e))}
spp213<-function(a,b,c,d,e) {rnorm(1, 10000, 4000)+(0*(a+b+c+d+e))}
spp214<-function(a,b,c,d,e) {rnorm(1, 10000, 400)+(0*(a+b+c+d+e))}
spp215<-function(a,b,c,d,e) {rnorm(1, 10000, 40)+(0*(a+b+c+d+e))}
spp216<-function(a,b,c,d,e) {rnorm(1, 10000, 4)+(0*(a+b+c+d+e))}
spp217<-function(a,b,c,d,e) {rnorm(1, 50000, 40000)+(0*(a+b+c+d+e))}
spp218<-function(a,b,c,d,e) {rnorm(1, 50000, 4000)+(0*(a+b+c+d+e))}
spp219<-function(a,b,c,d,e) {rnorm(1, 50000, 400)+(0*(a+b+c+d+e))}
spp220<-function(a,b,c,d,e) {rnorm(1, 50000, 40)+(0*(a+b+c+d+e))}
spp221<-function(a,b,c,d,e) {rnorm(1, 50000, 4)+(0*(a+b+c+d+e))}
spp222<-function(a,b,c,d,e) {rnorm(1, 1000000, 80000)+(0*(a+b+c+d+e))}
spp223<-function(a,b,c,d,e) {rnorm(1, 1000000, 8000)+(0*(a+b+c+d+e))}
spp224<-function(a,b,c,d,e) {rnorm(1, 1000000, 800)+(0*(a+b+c+d+e))}
spp225<-function(a,b,c,d,e) {rnorm(1, 1000000, 80)+(0*(a+b+c+d+e))}
spp226<-function(a,b,c,d,e) {rnorm(1, 1000000, 8)+(0*(a+b+c+d+e))}
spp227<-function(a,b,c,d,e) {rnorm(1, 60000, 40000)+(0*(a+b+c+d+e))}
spp228<-function(a,b,c,d,e) {rnorm(1, 50000, 40000)+(0*(a+b+c+d+e))}
spp229<-function(a,b,c,d,e) {rnorm(1, 40000, 40000)+(0*(a+b+c+d+e))}
spp230<-function(a,b,c,d,e) {rnorm(1, 30000, 40000)+(0*(a+b+c+d+e))}
spp231<-function(a,b,c,d,e) {rnorm(1, 20000, 40000)+(0*(a+b+c+d+e))}

# medium

spp232<-function(a,b,c,d,e) {rnorm(1, 1000, 10000)+(0*(a+b+c+d+e))}
spp233<-function(a,b,c,d,e) {rnorm(1, 1000, 1000)+(0*(a+b+c+d+e))}
spp234<-function(a,b,c,d,e) {rnorm(1, 1000, 100)+(0*(a+b+c+d+e))}
spp235<-function(a,b,c,d,e) {rnorm(1, 1000, 10)+(0*(a+b+c+d+e))}
spp236<-function(a,b,c,d,e) {rnorm(1, 1000, 1)+(0*(a+b+c+d+e))}
spp237<-function(a,b,c,d,e) {rnorm(1, 10000, 100000)+(0*(a+b+c+d+e))}
spp238<-function(a,b,c,d,e) {rnorm(1, 10000, 10000)+(0*(a+b+c+d+e))}
spp239<-function(a,b,c,d,e) {rnorm(1, 10000, 1000)+(0*(a+b+c+d+e))}
spp240<-function(a,b,c,d,e) {rnorm(1, 10000, 100)+(0*(a+b+c+d+e))}
spp241<-function(a,b,c,d,e) {rnorm(1, 10000, 10)+(0*(a+b+c+d+e))}
spp242<-function(a,b,c,d,e) {rnorm(1, 10000, 1)+(0*(a+b+c+d+e))}
spp243<-function(a,b,c,d,e) {rnorm(1, 100, 50)+(0*(a+b+c+d+e))}
spp244<-function(a,b,c,d,e) {rnorm(1, 200, 50)+(0*(a+b+c+d+e))}
spp245<-function(a,b,c,d,e) {rnorm(1, 300, 50)+(0*(a+b+c+d+e))}
spp246<-function(a,b,c,d,e) {rnorm(1, 400, 50)+(0*(a+b+c+d+e))}
spp247<-function(a,b,c,d,e) {rnorm(1, 500, 50)+(0*(a+b+c+d+e))}
spp248<-function(a,b,c,d,e) {rnorm(1, 600, 50)+(0*(a+b+c+d+e))}
spp249<-function(a,b,c,d,e) {rnorm(1, 700, 50)+(0*(a+b+c+d+e))}
spp250<-function(a,b,c,d,e) {rnorm(1, 800, 50)+(0*(a+b+c+d+e))}
spp251<-function(a,b,c,d,e) {rnorm(1, 900, 50)+(0*(a+b+c+d+e))}
spp252<-function(a,b,c,d,e) {rnorm(1, 100, 1000)+(0*(a+b+c+d+e))}
spp253<-function(a,b,c,d,e) {rnorm(1, 100, 100)+(0*(a+b+c+d+e))}
spp254<-function(a,b,c,d,e) {rnorm(1, 100, 10)+(0*(a+b+c+d+e))}
spp255<-function(a,b,c,d,e) {rnorm(1, 100, 1)+(0*(a+b+c+d+e))}
spp256<-function(a,b,c,d,e) {rnorm(1, 100, 100)+(0*(a+b+c+d+e))}
spp257<-function(a,b,c,d,e) {rnorm(1, 200, 200)+(0*(a+b+c+d+e))}
spp258<-function(a,b,c,d,e) {rnorm(1, 300, 300)+(0*(a+b+c+d+e))}
spp259<-function(a,b,c,d,e) {rnorm(1, 600, 600)+(0*(a+b+c+d+e))}
spp260<-function(a,b,c,d,e) {rnorm(1, 700, 700)+(0*(a+b+c+d+e))}
spp261<-function(a,b,c,d,e) {rnorm(1, 800, 800)+(0*(a+b+c+d+e))}
spp262<-function(a,b,c,d,e) {rnorm(1, 900, 900)+(0*(a+b+c+d+e))}

# Low

spp263<-function(a,b,c,d,e) {rnorm(100, 1000, 5)+(0*(a+b+c+d+e))}
spp264<-function(a,b,c,d,e) {rnorm(100, 100, 5)+(0*(a+b+c+d+e))}
spp265<-function(a,b,c,d,e) {rnorm(100, 10, 5)+(0*(a+b+c+d+e))}
spp266<-function(a,b,c,d,e) {rnorm(100, 1, 5)+(0*(a+b+c+d+e))}
spp267<-function(a,b,c,d,e) {rnorm(10, 1000, 5)+(0*(a+b+c+d+e))}
spp268<-function(a,b,c,d,e) {rnorm(10, 100, 5)+(0*(a+b+c+d+e))}
spp269<-function(a,b,c,d,e) {rnorm(10, 10, 5)+(0*(a+b+c+d+e))}
spp270<-function(a,b,c,d,e) {rnorm(10, 1, 5)+(0*(a+b+c+d+e))}
spp271<-function(a,b,c,d,e) {rnorm(1, 1000, 5)+(0*(a+b+c+d+e))}
spp272<-function(a,b,c,d,e) {rnorm(1, 100, 5)+(0*(a+b+c+d+e))}
spp273<-function(a,b,c,d,e) {rnorm(1, 10, 5)+(0*(a+b+c+d+e))}
spp274<-function(a,b,c,d,e) {rnorm(1, 1, 5)+(0*(a+b+c+d+e))}
spp275<-function(a,b,c,d,e) {rnorm(2, 1, 5)+(0*(a+b+c+d+e))}
spp276<-function(a,b,c,d,e) {rnorm(3, 1, 5)+(0*(a+b+c+d+e))}
spp277<-function(a,b,c,d,e) {rnorm(4, 1, 5)+(0*(a+b+c+d+e))}
spp278<-function(a,b,c,d,e) {rnorm(5, 1, 5)+(0*(a+b+c+d+e))}
spp279<-function(a,b,c,d,e) {rnorm(6, 1, 5)+(0*(a+b+c+d+e))}
spp280<-function(a,b,c,d,e) {rnorm(7, 1, 5)+(0*(a+b+c+d+e))}
spp281<-function(a,b,c,d,e) {rnorm(8, 1, 5)+(0*(a+b+c+d+e))}
spp282<-function(a,b,c,d,e) {rnorm(9, 1, 5)+(0*(a+b+c+d+e))}
spp283<-function(a,b,c,d,e) {rnorm(1, 1, 5)+(0*(a+b+c+d+e))}
spp284<-function(a,b,c,d,e) {rnorm(2, 2, 5)+(0*(a+b+c+d+e))}
spp285<-function(a,b,c,d,e) {rnorm(3, 3, 5)+(0*(a+b+c+d+e))}
spp286<-function(a,b,c,d,e) {rnorm(4, 4, 5)+(0*(a+b+c+d+e))}
spp287<-function(a,b,c,d,e) {rnorm(5, 5, 5)+(0*(a+b+c+d+e))}
spp288<-function(a,b,c,d,e) {rnorm(6, 6, 5)+(0*(a+b+c+d+e))}
spp289<-function(a,b,c,d,e) {rnorm(7, 7, 5)+(0*(a+b+c+d+e))}
spp290<-function(a,b,c,d,e) {rnorm(8, 8, 5)+(0*(a+b+c+d+e))}
spp291<-function(a,b,c,d,e) {rnorm(9, 9, 5)+(0*(a+b+c+d+e))}
spp292<-function(a,b,c,d,e) {rnorm(10, 10, 5)+(0*(a+b+c+d+e))}
spp293<-function(a,b,c,d,e) {rnorm(20, 20, 5)+(0*(a+b+c+d+e))}
spp294<-function(a,b,c,d,e) {rnorm(30, 30, 5)+(0*(a+b+c+d+e))}
spp295<-function(a,b,c,d,e) {rnorm(40, 40, 5)+(0*(a+b+c+d+e))}
spp296<-function(a,b,c,d,e) {rnorm(50, 50, 5)+(0*(a+b+c+d+e))}
spp297<-function(a,b,c,d,e) {rnorm(60, 60, 5)+(0*(a+b+c+d+e))}
spp298<-function(a,b,c,d,e) {rnorm(70, 70, 5)+(0*(a+b+c+d+e))}
spp299<-function(a,b,c,d,e) {rnorm(80, 80, 5)+(0*(a+b+c+d+e))}
spp300<-function(a,b,c,d,e) {rnorm(90, 90, 5)+(0*(a+b+c+d+e))}
