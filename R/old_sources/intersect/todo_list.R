todoListInput <- function(inputId, leftLabel, rightLabel, leftChoices, rightChoices,
                          size = 5, multiple = FALSE) {
  
  leftChoices <- lapply(leftChoices, tags$option)
  rightChoices <- lapply(rightChoices, tags$option)
  
  if (multiple)
    multiple <- "multiple"
  else
    multiple <- NULL
  
  
  cutsom_html = HTML(readLines("www/todo_list.html"))
  cutsom_html = gsub("SUBID", inputId, cutsom_html)
  cutsom_html = gsub("SUBLABEL", leftLabel, cutsom_html)
  
  outhtml = tagList(
    singleton(tags$head(
      tags$script(src="todo_list.js"),
      tags$style(type="text/css",
                 HTML('<link rel="stylesheet" href="todo_list.css">')
      )
    )),
    cutsom_html
  )
  
  return(outhtml)
}


'<head>
  <script src="todo_list.js"></script>
    </head>
    <div id="main">
      <h3>Tasks</h3>
      <div class="container" data-bind="sortable: tasks">
        <div class="item">
          <span data-bind="visible: !$root.isTaskSelected($data)">
            <a href="#" data-bind="text: name, click: $root.selectedTask"></a>
              </span>
              <span data-bind="visibleAndSelect: $root.isTaskSelected($data)">
                <input data-bind="value: name, event: { blur: $root.clearTask }" />
                  </span>  
                  </div>
                  </div>
                  <a href="#" data-bind="click: addTask">Add Task</a>
                    </div>
                    
                    <div id="results">
                      <h3>Tasks</h3>
                      <ul data-bind="foreach: tasks">
                        <li data-bind="text: name"></li>
                          </ul>
                          </div>'

registerInputHandler("shinyjsexamples.chooser", function(data, ...) {
  if (is.null(data))
    NULL
  else
    list(left=as.character(data$left), right=as.character(data$right))
}, force = TRUE)
