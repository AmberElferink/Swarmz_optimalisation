dear imgui
=====
[![Build Status](https://api.travis-ci.com/ocornut/imgui.svg?branch=master)](https://travis-ci.com/ocornut/imgui)
[![Coverity Status](https://scan.coverity.com/projects/4720/badge.svg)](https://scan.coverity.com/projects/4720)

<sub>(This library is available under a free and permissive license, but needs financial support to sustain its continued improvements. In addition to maintenance and stability there are many desirable features yet to be added. If your company is using dear imgui, please consider reaching out. If you are an individual using dear imgui, please consider supporting the project via Patreon or PayPal.)</sub>

Businesses: support continued development via invoiced technical support, maintenance, sponsoring contracts:
<br>&nbsp;&nbsp;_E-mail: contact @ dearimgui dot org_

Individuals/hobbyists: support continued maintenance and development via the monthly Patreon:
<br>&nbsp;&nbsp;[![Patreon](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/patreon_02.png)](http://www.patreon.com/imgui)

Individuals/hobbyists: support continued maintenance and development via PayPal:
<br>&nbsp;&nbsp;[![PayPal](https://www.paypalobjects.com/en_US/i/btn/btn_donate_LG.gif)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=WGHNC6MBFLZ2S)

----

Dear ImGui is a **bloat-free graphical user interface library for C++**. It outputs optimized vertex buffers that you can render anytime in your 3D-pipeline enabled application. It is fast, portable, renderer agnostic and self-contained (no external dependencies).

Dear ImGui is designed to **enable fast iterations** and to **empower programmers** to create **content creation tools and visualization / debug tools** (as opposed to UI for the average end-user). It favors simplicity and productivity toward this goal, and lacks certain features normally found in more high-level libraries.

Dear ImGui is particularly suited to integration in games engine (for tooling), real-time 3D applications, fullscreen applications, embedded applications, or any applications on consoles platforms where operating system features are non-standard.

| [Usage](#usage) - [How it works](#how-it-works) - [Demo](#demo) - [Integration](#integration) |
:----------------------------------------------------------: |
| [Upcoming changes](#upcoming-changes) - [Gallery](#gallery) - [Support, FAQ](#support-frequently-asked-questions-faq) -  [How to help](#how-to-help) - [Sponsors](#sponsors) - [Credits](#credits) - [License](#license) |
| [Wiki](https://github.com/ocornut/imgui/wiki) - [Language & frameworks bindings](https://github.com/ocornut/imgui/wiki/Bindings) - [Software using Dear ImGui](https://github.com/ocornut/imgui/wiki/Software-using-dear-imgui) - [User quotes](https://github.com/ocornut/imgui/wiki/Quotes) |

### Usage

Dear ImGui is self-contained within a few files that you can easily copy and compile into your application/engine:
- imgui.cpp
- imgui.h
- imgui_demo.cpp
- imgui_draw.cpp
- imgui_widgets.cpp
- imgui_internal.h
- imconfig.h (empty by default, user-editable)
- imstb_rectpack.h
- imstb_textedit.h
- imstb_truetype.h

No specific build process is required. You can add the .cpp files to your project or #include them from an existing file.

Backends for a variety of graphics api and rendering platforms along with example applications are provided in the [examples/](https://github.com/ocornut/imgui/tree/master/examples) folder.

The backend passes mouse/keyboard/gamepad inputs and variety of settings to Dear ImGui, and is in charge of rendering the resulting vertices. After Dear ImGui is setup in your application, you can use it from \_anywhere\_ in your program loop:

Code:
```cp
ImGui::Text("Hello, world %d", 123);
if (ImGui::Button("Save"))
    MySaveFunction();
ImGui::InputText("string", buf, IM_ARRAYSIZE(buf));
ImGui::SliderFloat("float", &f, 0.0f, 1.0f);
```
Result:
<br>![sample code output](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/v160/code_sample_02.png)
<br>_(settings: Dark style (left), Light style (right) / Font: Roboto-Medium, 16px / Rounding: 5)_

Code:
```cpp
// Create a window called "My First Tool", with a menu bar.
ImGui::Begin("My First Tool", &my_tool_active, ImGuiWindowFlags_MenuBar);
if (ImGui::BeginMenuBar())
{
    if (ImGui::BeginMenu("File"))
    {
        if (ImGui::MenuItem("Open..", "Ctrl+O")) { /* Do stuff */ }
        if (ImGui::MenuItem("Save", "Ctrl+S"))   { /* Do stuff */ }
        if (ImGui::MenuItem("Close", "Ctrl+W"))  { my_tool_active = false; }
        ImGui::EndMenu();
    }
    ImGui::EndMenuBar();
}

// Edit a color (stored as ~4 floats)
ImGui::ColorEdit4("Color", my_color);

// Plot some values
const float my_values[] = { 0.2f, 0.1f, 1.0f, 0.5f, 0.9f, 2.2f };
ImGui::PlotLines("Frame Times", my_values, IM_ARRAYSIZE(my_values));

// Display contents in a scrolling region
ImGui::TextColored(ImVec4(1,1,0,1), "Important Stuff");
ImGui::BeginChild("Scrolling");
for (int n = 0; n < 50; n++)
    ImGui::Text("%04d: Some text", n);
ImGui::EndChild();
ImGui::End();
```
Result:
<br>![sample code output](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/v160/code_sample_03_color.gif)

Dear ImGui allows you **create elaborate tools** as well as very short-lived ones. On the extreme side of short-liveness: using the Edit&Continue (hot code reload) feature of modern compilers you can add a few widgets to tweaks variables while your application is running, and remove the code a minute later! Dear ImGui is not just for tweaking values. You can use it to trace a running algorithm by just emitting text commands. You can use it along with your own reflection data to browse your dataset live. You can use it to expose the internals of a subsystem in your engine, to create a logger, an inspection tool, a profiler, a debugger, an entire game making editor/framework, etc.

### How it works

Check out the Wiki's [About the IMGUI paradigm](https://github.com/ocornut/imgui/wiki#About-the-IMGUI-paradigm) section if you want to understand the core principles behind the IMGUI paradigm. An IMGUI tries to minimize superfluous state duplication, state synchronization and state retention from the user's point of view. It is less error prone (less code and less bugs) than traditional retained-mode interfaces, and lends itself to create dynamic user interfaces.

Dear ImGui outputs vertex buffers and command lists that you can easily render in your application. The number of draw calls and state changes required to render them is fairly small. Because Dear ImGui doesn't know or touch graphics state directly, you can call its functions  anywhere in your code (e.g. in the middle of a running algorithm, or in the middle of your own rendering process). Refer to the sample applications in the examples/ folder for instructions on how to integrate dear imgui with your existing codebase.

_A common misunderstanding is to mistake immediate mode gui for immediate mode rendering, which usually implies hammering your driver/GPU with a bunch of inefficient draw calls and state changes as the gui functions are called. This is NOT what Dear ImGui does. Dear ImGui outputs vertex buffers and a small list of draw calls batches. It never touches your GPU directly. The draw call batches are decently optimal and you can render them later, in your app or even remotely._

### Demo

Calling the `ImGui::ShowDemoWindow()` function will create a demo window showcasing variety of features and examples. The code is always available for reference in `imgui_demo.cpp`.

![screenshot demo](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/v167/v167-misc.png)

You should be able to build the examples from sources (tested on Windows/Mac/Linux). If you don't, let me know! If you want to have a quick look at some Dear ImGui features, you can download Windows binaries of the demo app here:
- [imgui-demo-binaries-20190715.zip](http://www.dearimgui.org/binaries/imgui-demo-binaries-20190715.zip) (Windows binaries, 1.72 WIP, built 2019/07/15, master branch, 5 executables)

The demo applications are not DPI aware so expect some blurriness on a 4K screen. For DPI awareness in your application, you can load/reload your font at different scale, and scale your style with `style.ScaleAllSizes()`.

### Integration

On most platforms and when using C++, **you should be able to use a combination of the [imgui_impl_xxxx](https://github.com/ocornut/imgui/tree/master/examples) files without modification** (e.g. `imgui_impl_win32.cpp` + `imgui_impl_dx11.cpp`). If your engine supports multiple platforms, consider using more of the imgui_impl_xxxx files instead of rewriting them: this will be less work for you and you can get Dear ImGui running immediately. You can _later_ decide to rewrite a custom binding using your custom engine functions if you wish so.

Integrating Dear ImGui within your custom engine is a matter of 1) wiring mouse/keyboard/gamepad inputs 2) uploading one texture to your GPU/render engine 3) providing a render function that can bind textures and render textured triangles. The [examples/](https://github.com/ocornut/imgui/tree/master/examples) folder is populated with applications doing just that. If you are an experienced programmer at ease with those concepts, it should take you less than two hours to integrate Dear ImGui in your custom engine. **Make sure to spend time reading the FAQ, comments, and one of the examples/ application!**

Officially maintained bindings (in repository):
- Renderers: DirectX9, DirectX10, DirectX11, DirectX12, OpenGL (legacy), OpenGL3/ES/ES2 (modern), Vulkan, Metal.
- Platforms: GLFW, SDL2, Win32, Glut, OSX.
- Frameworks: Emscripten, Allegro5, Marmalade.

Third-party bindings (see [Bindings](https://github.com/ocornut/imgui/wiki/Bindings/) page):
- Languages: C, C#/.Net, ChaiScript, D, Go, Haxe/hxcpp, Java, JavaScript, Julia, Lua, Odin, Pascal, PureBasic, Python, Ruby, Rust, Swift...
- Frameworks: bsf, Cinder, Cocos2d-x, Diligent Engine, Flexium, GML/GameMakerStudio2, Irrlicht, Ogre, OpenFrameworks, OpenSceneGraph/OSG, ORX, px_render, L??VE+Lua, Magnum, NanoRT, Qt, QtDirect3D, SFML, Software Rasterizers, Unreal Engine 4...
- Note that C bindings ([cimgui](https://github.com/cimgui/cimgui)) are auto-generated, you can use its json/lua output to generate bindings for other languages.

Also see [Wiki](https://github.com/ocornut/imgui/wiki) for more links and ideas.

### Upcoming Changes

Some of the goals for 2019 are:
- Finish work on docking, tabs. (see [#2109](https://github.com/ocornut/imgui/issues/2109), in public [docking](https://github.com/ocornut/imgui/tree/docking) branch looking for feedback)
- Finish work on multiple viewports / multiple OS windows. (see [#1542](https://github.com/ocornut/imgui/issues/1542), in public [docking](https://github.com/ocornut/imgui/tree/docking) branch looking for feedback)
- Finish work on gamepad/keyboard controls. (see [#787](https://github.com/ocornut/imgui/issues/787))
- Add an automation and testing system, both to test the library and end-user apps. (see [#435](https://github.com/ocornut/imgui/issues/435))
- Make Columns better. They are currently pretty terrible! New Tables API coming Q4 2019!
- Make the examples look better, improve styles, improve font support, make the examples hi-DPI and multi-DPI aware.

### Gallery

For more user-submitted screenshots of projects using Dear ImGui, check out the [Gallery Threads](https://github.com/ocornut/imgui/issues/2529)!

Custom engine
[![screenshot game](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/v149/gallery_TheDragonsTrap-01-thumb.jpg)](https://cloud.githubusercontent.com/assets/8225057/20628927/33e14cac-b329-11e6-80f6-9524e93b048a.png)

Custom engine
[![screenshot tool](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/v160/editor_white_preview.jpg)](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/v160/editor_white.png)

[Tracy Profiler](https://bitbucket.org/wolfpld/tracy)
![tracy profiler](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/v167/tracy_profiler.png)

### Wiki

See [Wiki](https://github.com/ocornut/imgui/wiki) for many links, references, articles.

See [Articles about the IMGUI paradigm](https://github.com/ocornut/imgui/wiki#Articles-about-the-IMGUI-paradigm) to read/watch about the Immediate Mode GUI paradigm.

### Support, Frequently Asked Questions (FAQ)

If you are new to Dear ImGui and have issues with: compiling, linking, adding fonts, wiring inputs, running or displaying Dear ImGui: you can use [Discord server](https://discord.gg/NgJ4SEP) or [Discourse forums](https://discourse.dearimgui.org).

Otherwise, for any other questions, bug reports, requests, feedback, you may post on https://github.com/ocornut/imgui/issues. Please read and fill the New Issue template carefully.

Paid private support is available for business customers (E-mail: _contact @ dearimgui dot org_).

**Where is the documentation?**

 This library is poorly documented at the moment and expects of the user to be acquainted with C/C++.
 - Run the examples/ applications and explore them.
 - See demo code in imgui_demo.cpp and particularly the ImGui::ShowDemoWindow() function.
 - The demo covers most features of Dear ImGui, so you can read the code and see its output.
 - See documentation and comments at the top of imgui.cpp + effectively imgui.h.
 - Dozens of standalone example applications using e.g. OpenGL/DirectX are provided in the examples/ folder to explain how to integrate Dear ImGui with your own engine/application.
 - Your programming IDE is your friend, find the type or function declaration to find comments associated with it.
 - We obviously need better documentation! Consider contributing or becoming a [Patron](http://www.patreon.com/imgui) to promote this effort.

**Which version should I get?**

I occasionally tag [Releases](https://github.com/ocornut/imgui/releases) but it is generally safe and recommended to sync to master/latest. The library is fairly stable and regressions tend to be fixed fast when reported.

You may also peak at the [Multi-Viewport](https://github.com/ocornut/imgui/issues/1542) and [Docking](https://github.com/ocornut/imgui/issues/2109) features in the `docking` branch. Many projects are using this branch and it is kept in sync with master regularly.

**Who uses Dear ImGui?**

See the [Quotes](https://github.com/ocornut/imgui/wiki/Quotes) and [Software using dear imgui](https://github.com/ocornut/imgui/wiki/Software-using-dear-imgui) Wiki pages for a list of games/software which are publicly known to use dear imgui. Please add yours if you can!

**Why the odd dual naming, "Dear ImGui" vs "ImGui"?**

The library started its life as "ImGui" due to the fact that I didn't give it a proper name when I released 1.0 and had no particular expectation that it would take off. However, the term IMGUI (immediate-mode graphical user interface) was coined before and is being used in variety of other situations (e.g. Unity uses its own implementation of the IMGUI paradigm). To reduce this ambiguity without affecting existing codebases, I have decided on an alternate, longer name "Dear ImGui" that people can use to refer to this specific library. Please try to refer to this library as "Dear ImGui".

**How can I tell whether to dispatch mouse/keyboard to Dear ImGui or to my application?**
<br>**How can I display an image? What is ImTextureID, how does it works?** ([examples](https://github.com/ocornut/imgui/wiki/Image-Loading-and-Displaying-Examples))
<br>**Why are multiple widgets reacting when I interact with a single one? How can I have multiple widgets with the same label or with an empty label? A primer on labels and the ID Stack...**
<br>**How can I use my own math types instead of ImVec2/ImVec4?**
<br>**How can I load a different font than the default?**
<br>**How can I easily use icons in my application?**
<br>**How can I load multiple fonts?**
<br>**How can I display and input non-latin characters such as Chinese, Japanese, Korean, Cyrillic?** ([examples](https://github.com/ocornut/imgui/wiki/Loading-Font-Example))
<br>**How can I interact with standard C++ types (such as std::string and std::vector)?**
<br>**How can I use the drawing facilities without a Dear ImGui window? (using ImDrawList API)**
<br>**How can I use this without a mouse, without a keyboard or without a screen? (gamepad, input share, remote display)**
<br>**I integrated Dear ImGui in my engine and the text or lines are blurry..**
<br>**I integrated Dear ImGui in my engine and some elements are disappearing when I move windows around..**
<br>**How can I help?**

See the FAQ in [imgui.cpp](https://github.com/ocornut/imgui/blob/master/imgui.cpp) for answers.

**Can you create elaborate/serious tools with Dear ImGui?**

Yes. People have written game editors, data browsers, debuggers, profilers and all sort of non-trivial tools with the library. In my experience the simplicity of the API is very empowering. Your UI runs close to your live data. Make the tools always-on and everybody in the team will be inclined to create new tools (as opposed to more "offline" UI toolkits where only a fraction of your team effectively creates tools). The list of sponsors below is also an indicator that serious game teams have been using the library.

Dear ImGui is very programmer centric and the immediate-mode GUI paradigm might require you to readjust some habits before you can realize its full potential. Dear ImGui is about making things that are simple, efficient and powerful.

**Can you reskin the look of Dear ImGui?**

You can alter the look of the interface to some degree: changing colors, sizes, padding, rounding, fonts. However, as Dear ImGui is designed and optimized to create debug tools, the amount of skinning you can apply is limited. There is only so much you can stray away from the default look and feel of the interface. Below is a screenshot from [LumixEngine](https://github.com/nem0/LumixEngine) with custom colors + a docking/tabs extension (both of which you can find in the Issues section and will eventually be merged):

![LumixEngine](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/v151/lumix-201710-rearranged.png)

**Why using C++ (as opposed to C)?**

Dear ImGui takes advantage of a few C++ languages features for convenience but nothing anywhere Boost insanity/quagmire. Dear ImGui does NOT require C++11 so it can be used with most old C++ compilers. Dear ImGui doesn't use any C++ header file. Language-wise, function overloading and default parameters are used to make the API easier to use and code more terse. Doing so I believe the API is sitting on a sweet spot and giving up on those features would make the API more cumbersome. Other features such as namespace, constructors and templates (in the case of the ImVector<> class) are also relied on as a convenience.

There is an auto-generated [c-api for Dear ImGui (cimgui)](https://github.com/cimgui/cimgui) by Sonoro1234 and Stephan Dilly. It is designed for creating binding to other languages. If possible, I would suggest using your target language functionalities to try replicating the function overloading and default parameters used in C++ else the API may be harder to use. Also see [Bindings](https://github.com/ocornut/imgui/wiki/Bindings) for various third-party bindings.

How to help
-----------

**How can I help?**

- You may participate in the [Discord server](https://discord.gg/NgJ4SEP), [Discourse forums](https://discourse.dearimgui.org), GitHub [issues tracker](https://github.com/ocornut/imgui/issues).
- You may help with development and submit pull requests! Please understand that by submitting a PR you are also submitting a request for the maintainer to review your code and then take over its maintenance forever. PR should be crafted both in the interest in the end-users and also to ease the maintainer into understanding and accepting it.
- See [Help wanted](https://github.com/ocornut/imgui/wiki/Help-Wanted) on the [Wiki](https://github.com/ocornut/imgui/wiki/) for some more ideas.
- Have your company financially support this project.

**How can I help financing further development of Dear ImGui?**

Your contributions are keeping this project alive. The library is available under a free and permissive license, but continued maintenance and development are a full-time endeavor and I would like to grow the team. In addition to maintenance and stability there are many desirable features yet to be added. If your company is using dear imgui, please consider reaching out for invoiced technical support and maintenance contracts. If you are an individual using dear imgui, please consider supporting the project via Patreon or PayPal. Thank you!

Businesses: support continued development via invoiced technical support, maintenance, sponsoring contracts:
<br>&nbsp;&nbsp;_E-mail: contact @ dearimgui dot org_

Individuals/hobbyists: support continued maintenance and development via the monthly Patreon:
<br>&nbsp;&nbsp;[![Patreon](https://raw.githubusercontent.com/wiki/ocornut/imgui/web/patreon_02.png)](http://www.patreon.com/imgui)

Individuals/hobbyists: support continued maintenance and development via PayPal:
<br>&nbsp;&nbsp;[![PayPal](https://www.paypalobjects.com/en_US/i/btn/btn_donate_LG.gif)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=WGHNC6MBFLZ2S)

### Sponsors

Ongoing Dear ImGui development is financially supported by users and private sponsors, recently:

*Platinum-chocolate sponsors*
- Blizzard Entertainment
- Google

*Double-chocolate sponsors*
- Media Molecule, Mobigame, Aras Pranckevi??ius, Greggman, DotEmu, Nadeo, Supercell, Runner, Aiden Koss, Kylotonn.

*Salty caramel supporters*
- Remedy Entertainment, Recognition Robotics, ikrima, Geoffrey Evans, Mercury Labs, Singularity Demo Group, Lionel Landwerlin, Ron Gilbert, Brandon Townsend, G3DVu, Cort Stratton, drudru, Harfang 3D, Jeff Roberts, Rainway inc, Ondra Voves, Mesh Consultants, Unit 2 Games, Neil Bickford, Bill Six, Graham Manders.

*Caramel supporters*
- Jerome Lanquetot, Daniel Collin, Ctrl Alt Ninja, Neil Henning, Neil Blakey-Milner, Aleksei, NeiloGD, Eric, Game Atelier, Vincent Hamm, Morten Skaaning, Colin Riley, Sergio Gonzales, Andrew Berridge, Roy Eltham, Game Preservation Society, Josh Faust, Martin Donlon, Codecat, Doug McNabb, Emmanuel Julien, Guillaume Chereau, Jeffrey Slutter, Jeremiah Deckard, r-lyeh, Nekith, Joshua Fisher, Malte Hoffmann, Mustafa Karaalioglu, Merlyn Morgan-Graham, Per Vognsen, Fabian Giesen, Jan Staubach, Matt Hargett, John Shearer, Jesse Chounard, kingcoopa, Jonas Bernemann, Johan Andersson, Michael Labbe, Tomasz Golebiowski, Louis Schnellbach, Jimmy Andrews, Bojan Endrovski, Robin Berg Pettersen, Rachel Crawford, Andrew Johnson, Sean Hunter, Jordan Mellow, Nefarius Software Solutions, Laura Wieme, Robert Nix, Mick Honey, Steven Kah Hien Wong, Bartosz Bielecki, Oscar Penas, A M, Liam Moynihan, Artometa, Mark Lee, Dimitri Diakopoulos, Pete Goodwin, Johnathan Roatch, nyu lea, Oswald Hurlem, Semyon Smelyanskiy, Le Bach, Jeong MyeongSoo, Chris Matthews, Astrofra, Frederik De Bleser, Anticrisis, Matt Reyer.

And all other past and present supporters; THANK YOU!
(Please contact me if you would like to be added or removed from this list)

Credits
-------

Developed by [Omar Cornut](http://www.miracleworld.net) and every direct or indirect contributors to the GitHub. The early version of this library was developed with the support of [Media Molecule](http://www.mediamolecule.com) and first used internally on the game [Tearaway](http://tearaway.mediamolecule.com) (Vita).

I first discovered the IMGUI paradigm at [Q-Games](http://www.q-games.com) where Atman Binstock had dropped his own simple implementation in the codebase, which I spent quite some time improving and thinking about. It turned out that Atman was exposed to the concept directly by working with Casey. When I moved to Media Molecule I rewrote a new library trying to overcome the flaws and limitations of the first one I've worked with. It became this library and since then I have spent an unreasonable amount of time iterating and improving it.

Embeds [ProggyClean.ttf](http://upperbounds.net) font by Tristan Grimmer (MIT license).

Embeds [stb_textedit.h, stb_truetype.h, stb_rect_pack.h](https://github.com/nothings/stb/) by Sean Barrett (public domain).

Inspiration, feedback, and testing for early versions: Casey Muratori, Atman Binstock, Mikko Mononen, Emmanuel Briney, Stefan Kamoda, Anton Mikhailov, Matt Willis. And everybody posting feedback, questions and patches on the GitHub.

License
-------

Dear ImGui is licensed under the MIT License, see [LICENSE.txt](https://github.com/ocornut/imgui/blob/master/LICENSE.txt) for more information.
