# assets/

Drop a TrueType font here named `font.ttf`. The simulator's font loader
also probes a few common system locations automatically, so on most
machines no manual step is required.

Examples (the file itself is **gitignored** to avoid shipping a
license-encumbered binary):

```bash
# Linux
cp /usr/share/fonts/truetype/dejavu/DejaVuSans.ttf assets/font.ttf

# Windows (PowerShell)
Copy-Item C:\Windows\Fonts\segoeui.ttf assets\font.ttf

# macOS
cp /System/Library/Fonts/Supplemental/Arial.ttf assets/font.ttf
```
