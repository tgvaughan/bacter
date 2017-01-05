function sidebarUpdate() {
    var origVPos = $("#sidebar").parent().position().top;
    var origWidth = $("#sidebar").parent().width();

    if ($(this).scrollTop() > origVPos-20) {
        $("#sidebar").addClass("sidebar-fixed");
        $("#sidebar").width(origWidth);
    } else {
        $("#sidebar").removeClass("sidebar-fixed");
        $("#sidebar").width("");
    }
}
$(window).scroll(sidebarUpdate);
$(window).resize(sidebarUpdate);
