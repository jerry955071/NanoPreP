template="""\
<html>
    <head><meta charset="utf-8" /></head>
    <style type="text/css">
        .gap-20 { 
            width:100%; 
            height:20px; 
        }
        .gap-10 { 
            width:100%; 
            height:10px; 
        } 
        #title {
            font-family: Sans-serif;
            font-size: 24px;
            font-weight: 600;
        }
        #text {
            font-family: Sans-serif;
            font-size: 16px;
            font-weight: 600;
        }
        #comment {
            font-family: Sans-serif;
            font-size: 12px;
        }
        #code {
            font-family: Monospace;
            font-size: 14px;
        }
    </style>
    <body>
        <div id="title">NanoPreP-optimize</div>
        <div class="gap-20"></div>
        <div>%(plot)</div>
        <div id="text">Command line:</div>
        <div id="code">nanoprep-optimize %(args)</div>
        <div class="gap-10"></div>
        <div id="text">Parameters:</div>
        <div id="code">%(params)</div>
    </body>
</html>
"""