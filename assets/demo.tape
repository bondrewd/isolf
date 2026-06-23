# isolf terminal demo.
#
# Records the build output (gradient title + twin-panel summary) to assets/demo.gif.
# Needs `isolf` on your PATH and vhs (https://github.com/charmbracelet/vhs), then:
#
#     vhs assets/demo.tape

Output assets/demo.gif

Set Shell "bash"
Set FontFamily "JuliaMono"
Set FontSize 15
Set Width 1180
Set Height 820
Set Padding 24
Set Theme "Catppuccin Mocha"

Hide
Type "clear" Enter
Show

Type "isolf --upper 'DOPC=3,DPPC=1' --lower 'DOPC=1,DPPC=1' --membrane 18 --inp --vmd --out demo"
Sleep 700ms
Enter
Sleep 4s
